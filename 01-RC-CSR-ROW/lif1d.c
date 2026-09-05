#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <getopt.h>
#include <stdlib.h>
#include <sys/time.h>
#include <dpu.h>
#include <dpu_log.h>
#include <assert.h>
#include <unistd.h>

#include "./header_files/matrix.h"
#include "./header_files/common.h"
#include "./header_files/partition.h"
#include "./header_files/timer.h"

#ifndef DPU_BINARY
#define DPU_BINARY "./dpu"
#endif
#define NR_TASKELTS 16

#define DPU_CAPACITY (64 << 20)

#define DEF_DT			(1.0e-04)
#define DEF_MU			(1.0)
#define DEF_UTH			(0.98)
#define	DEF_S_MIN		(0.7)
#define DEF_S_MAX		(0.7)
#define DEF_SAMPLING		(100L)
#define DEF_SIM_TIME		(20.0)
#define DEF_TTRANSIENT		(-1L)

static struct option long_options[] =
{
 {"dt",        required_argument, 0, 'a'},
 {"transient", required_argument, 0, 'e'},
 {"time",      required_argument, 0, 'd'},
 {"s_min",     required_argument, 0, 'f'},
 {"s_max",     required_argument, 0, 'g'},
 {"sampling",  required_argument, 0, 'h'},
 {"spfile",    required_argument, 0, 'i'},
 {0, 0, 0, 0}
};

//added 1.
static struct partition_info_t *part_info;


	struct dpu_info_t {
    uint32_t rows_per_dpu;
    uint32_t rows_per_dpu_pad;
    uint32_t prev_rows_dpu;
    uint32_t nnz;
    uint32_t nnz_pad;
};
struct dpu_info_t *dpu_info;

int main(int argc, char *argv[])
{
	FILE		*output1, *output2; //outfiles after of the calculation
	long		n; //number of neurons
	long		i, j; //individual neurons
	long		it; //number of iterations
	long		sampling;
	double		dt; //time diference Δt
	long		time_steps_per_sec;
	double		sim_time;
	long		ttransient; //time to stabylize
	long		total_time_steps;
	long		*divide; //Number of neighours per neuron
	double		uth; //thresshold
	double		mu;
	double		s_min; //lower bound of values
	double		s_max; //upper bound of values
	double		*u, *u_next, *omega /*frequency of spike*/, *omega1 /*times of spike*/, *sum_s /* to athrisma ths eksiswshs*/, *temp; //u = u(t), u_next = u(t+Δt), temp = for u=u_next.
	double		time;
	struct timeval	global_start, global_end, IO_start, IO_end;
	double		global_usec, IO_usec = 0.0;
	int		c, option_index;
	char		*end_ptr, transa = 'N';
	char		*filename0 = NULL, filename1[NAME_MAX + 1], filename2[NAME_MAX + 1];
	static struct COOMatrix*	s_coo; // Attention COO matrix
	static struct CSRMatrix*	s; //Attention CSR matrix
	unsigned int k; //In christinas this variable has the name "i"
	unsigned int cur_rows = 0;
	double *y;
	double average_run_time = 0;
	double total_average_run_time = 0;
	double average_vector_allocation_time = 0;
	double total_vector_allocation_time = 0;
	double average_output_gather_time = 0;
	double total_output_gather_time=0;

	struct dpu_set_t dpu_set;
	struct dpu_set_t dpu;
	uint32_t nr_of_dpus;

	dt         = DEF_DT;
	mu         = DEF_MU;
	uth        = DEF_UTH;
	s_min      = DEF_S_MIN;
	s_max      = DEF_S_MAX;
	sampling   = DEF_SAMPLING;
	sim_time   = DEF_SIM_TIME;
	ttransient = DEF_TTRANSIENT;

	while (1) {
		c = getopt_long (argc, argv, "", long_options, &option_index);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 'a':
				dt = strtod(optarg, &end_ptr);
				if (*end_ptr != '\0') {
					printf("Option \"%s\": Invalid argument \"%s\".\n", long_options[option_index].name, optarg);
					exit(1);
				}
				if (dt <= 0.0) {
					printf("Option \"%s\": \"dt\" must be larger than zero.\n", long_options[option_index].name);
					exit(1);
				}
				break;
			case 'b':
				mu = strtod(optarg, &end_ptr);
				if (*end_ptr != '\0') {
					printf("Option \"%s\": Invalid argument \"%s\".\n", long_options[option_index].name, optarg);
					exit(1);
				}
				if (mu <= 0.0) {
					printf("Option \"%s\": \"mu\" must be larger than zero.\n", long_options[option_index].name);
					exit(1);
				}
				break;
			case 'c':
				uth = strtod(optarg, &end_ptr);
				if (*end_ptr != '\0') {
					printf("Option \"%s\": Invalid argument \"%s\".\n", long_options[option_index].name, optarg);
					exit(1);
				}
				if (uth <= 0.0) {
					printf("Option \"%s\": \"uth\" must be larger than zero.\n", long_options[option_index].name);
					exit(1);
				}
				break;
			case 'd':
				sim_time = strtod(optarg, &end_ptr);
				if (*end_ptr != '\0') {
					printf("Option \"%s\": Invalid argument \"%s\".\n", long_options[option_index].name, optarg);
					exit(1);
				}
				/*if (sim_time < 1) {
					printf("Option \"%s\": Total simulation time must be larger than zero.\n", long_options[option_index].name);
					exit(1);
				}*/
				break;
			case 'e':
				ttransient = strtol(optarg, &end_ptr, 10);
				if (*end_ptr != '\0') {
					printf("Option \"%s\": Invalid argument \"%s\".\n", long_options[option_index].name, optarg);
					exit(1);
				}
				if (ttransient < 0) {
					printf("Option \"%s\": \"ttransient\" must be larger or equal than zero.\n", long_options[option_index].name);
					exit(1);
				}
				break;
			case 'f':
				s_min = strtod(optarg, &end_ptr);
				if (*end_ptr != '\0') {
					printf("Option \"%s\": Invalid argument \"%s\".\n", long_options[option_index].name, optarg);
					exit(1);
				}
				break;
			case 'g':
				s_max = strtod(optarg, &end_ptr);
				if (*end_ptr != '\0') {
					printf("Option \"%s\": Invalid argument \"%s\".\n", long_options[option_index].name, optarg);
					exit(1);
				}
				break;
			case 'h':
				sampling = strtol(optarg, &end_ptr, 10);
				if (*end_ptr != '\0') {
					printf("Option \"%s\": Invalid argument \"%s\".\n", long_options[option_index].name, optarg);
					exit(1);
				}
				if (sampling < 1) {
					printf("Option \"%s\": \"sampling\" must be larger or equal than one.\n", long_options[option_index].name);
					exit(1);
				}
				break;
			case 'i':
				filename0 = optarg;
				break;
			case '?':
			default:
				exit(1);
				break;
		}
	}

	if (optind != argc) {
		printf("Unknown option \"%s\".\n", argv[optind]);
		exit(1);
	}

	if (filename0 == NULL) {
		printf("Please provide an input file for the connectivity matrix.\n");
		exit(1);
	}

	if (s_min > s_max) {
		printf("s_min (%17.15f) must be smaller or equal than s_max (%17.15f).\n", s_min, s_max);
		exit(1);
	}

	time_steps_per_sec = (long)(1.0 / dt);
	if (ttransient == DEF_TTRANSIENT) {
		ttransient = (sim_time * time_steps_per_sec) / 2;
	} else {
		ttransient *= time_steps_per_sec;
	}
	total_time_steps = sim_time * time_steps_per_sec;

	// Allocate DPUs and load binary
    DPU_ASSERT(dpu_alloc(NR_DPUS, NULL, &dpu_set));
    DPU_ASSERT(dpu_load(dpu_set, DPU_BINARY, NULL));
    DPU_ASSERT(dpu_get_nr_dpus(dpu_set, &nr_of_dpus));
    printf("[INFO] Allocated %d DPU(s)\n", nr_of_dpus);
    printf("[INFO] Allocated %d TASKLET(s) per DPU\n", NR_TASKELTS);

	s_coo = readCOOMatrix(filename0, s_min, s_max); 
	s = coo2csr(s_coo); 
	freeCOOMatrix(s_coo); 

	if (s->nrows != s->ncols) { 
		printf("Connectivity matrix must be square.\n");
		exit(1);
	}
	n = s->nrows; //n = number of neurons

	// Initialize partition data
    part_info = partition_init(nr_of_dpus, NR_TASKELTS);

    // Load-balance nnz across DPUs
    partition_by_row(s, part_info, nr_of_dpus);

	printf("Running simulation with following parameters:\n");
	printf("  Number of neurons   : %ld\n", n);
	printf("  Connectivity matrix : %s\n", filename0);
	printf("  Simulation time     : %f seconds (%ld time steps)\n", sim_time, total_time_steps);
	printf("  Transient time      : %ld seconds (%ld time steps)\n", ttransient / time_steps_per_sec, ttransient);
	printf("  Sampling rate       : %ld time steps\n", sampling);
	printf("  dt                  : %.1e seconds \n", dt);
	printf("  mu                  : %17.15f\n", mu);
	printf("  uth                 : %17.15f\n", uth);
	printf("  s_min               : %17.15f\n", s_min);
	printf("  s_max               : %17.15f\n", s_max);
	//name of first output file
	i = snprintf(filename1, NAME_MAX + 1, "spacetime_%07ld_%+6.4f_%+6.4f_%8.6f_%05f_%05ld_%8.6f_%8.6f_%07ld.out",
						n, s_min, s_max, dt, sim_time, ttransient / time_steps_per_sec, mu, uth, sampling);
	if (i >= NAME_MAX + 1) {
		printf("Filename to store space-time information is too long.\n");
		exit(1);
	}
	//name of second output file
	i = snprintf(filename2, NAME_MAX + 1, "omega_%07ld_%+6.4f_%+6.4f_%8.6f_%05f_%05ld_%8.6f_%8.6f_%07ld.out",
						n, s_min, s_max, dt, sim_time, ttransient / time_steps_per_sec, mu, uth, sampling);
	if (i >= NAME_MAX + 1) {
		printf("Filename to store omega information is too long.\n");
		exit(1);
	}
	//opening/creating files
	output1 = fopen(filename1, "w");
	if (output1 == NULL) {
		printf("Could not open file \"%s\".\n", filename1);
		exit(1);
	}

	output2 = fopen(filename2, "w");
	if (output2 == NULL) {
		printf("Could not open file \"%s\".\n", filename2);
		exit(1);
	}
	
	//allocate memeory for vector
	u = (double *)calloc(n, sizeof(double));
	if (u == NULL) {
		printf("Could not allocate memory for \"u\".\n");
		exit(1);
	}
	
	u_next = (double *)calloc(n, sizeof(double));
	if (u_next == NULL) {
		printf("Could not allocate memory for \"u_next\".\n");
		exit(1);
	}

	omega = (double *)calloc(n, sizeof(double));
	if (omega == NULL) {
		printf("Could not allocate memory for \"omega\".\n");
		exit(1);
	}

	omega1 = (double *)calloc(n, sizeof(double));
	if (omega1 == NULL) {
		printf("Could not allocate memory for \"omega1\".\n");
		exit(1);
	}

	sum_s = (double *)calloc(n, sizeof(double));
	if (sum_s == NULL) {
		printf("Could not allocate memory for \"sum_s\".\n");
		exit(1);
	}

	divide = (long *)calloc(n, sizeof(long));
	if (divide == NULL) {
		printf("Could not allocate memory for \"divide\".\n");
		exit(1);
	}

	/*
	 * Initialize elements of array u[i] with random numbers in the range [0, uth).
	 */
	for (i = 0; i < n; i++ ) {
		u[i] = uth * drand48();
	}
	
	// Initialize help data
    dpu_info = (struct dpu_info_t *) malloc(nr_of_dpus * sizeof(struct dpu_info_t)); 
    dpu_arguments_t *input_args = (dpu_arguments_t *) malloc(nr_of_dpus * sizeof(dpu_arguments_t));
    // Max limits for parallel transfers
    uint64_t max_rows_per_dpu = 0;
    uint64_t max_nnz_ind_per_dpu = 0;
    uint64_t max_nnz_val_per_dpu = 0;
    uint64_t max_rows_per_tasklet = 0;

#if 0
	printCSRMatrix(s);
#endif
//Before padding becasue reallocating s.
	for (i = 0; i < s->nrows; i++) {
		for (j = s->rowptr[i]; j < s->rowptr[i + 1]; j++) {
			sum_s[i] += s->values[j];
			/*
			 * If there is a connection between neuron i and neuron j
			 * update the number of neighbours for neuron i.
			 */
			if ((i != s->colind[j]) && (s->values[j] != 0.0)) {
				divide[i]++;
			}
		}
		if (divide[i] == 0) {
			divide[i]++;
		}
	}

    // Timer for measurements
    Timer timer;

    k = 0;
    DPU_FOREACH(dpu_set, dpu, k) {
        // Find padding for rows and non-zero elements needed for CPU-DPU transfers
        uint32_t rows_per_dpu = part_info->row_split[k+1] - part_info->row_split[k];
        uint32_t prev_rows_dpu = part_info->row_split[k];

        // Pad data to be transfered for rows
        uint32_t rows_per_dpu_pad = rows_per_dpu + 1;
        if (rows_per_dpu_pad % (8 / byte_dt) != 0)
            rows_per_dpu_pad += ((8 / byte_dt) - (rows_per_dpu_pad % (8 / byte_dt)));
#if INT64 || FP64
        if (rows_per_dpu_pad % 2 == 1)
            rows_per_dpu_pad++;
#endif
        if (rows_per_dpu_pad > max_rows_per_dpu)
            max_rows_per_dpu = rows_per_dpu_pad;

        // Pad data to be transfered for nnzs
        unsigned int nnz, nnz_ind_pad, nnz_val_pad;
        nnz = s->rowptr[rows_per_dpu + prev_rows_dpu] - s->rowptr[prev_rows_dpu];
        if (nnz % 2 != 0)
            nnz_ind_pad = nnz + 1;
        else
            nnz_ind_pad = nnz;
        if (nnz % (8 / byte_dt) != 0)
            nnz_val_pad = nnz + ((8 / byte_dt) - (nnz % (8 / byte_dt)));
        else
            nnz_val_pad = nnz;

#if INT64 || FP64
        if (nnz_ind_pad % 2 == 1)
            nnz_ind_pad++;
        if (nnz_val_pad % 2 == 1)
            nnz_val_pad++;
#endif
        if (nnz_ind_pad > max_nnz_ind_per_dpu)
            max_nnz_ind_per_dpu = nnz_ind_pad;
        if (nnz_val_pad > max_nnz_val_per_dpu)
            max_nnz_val_per_dpu = nnz_val_pad;

        // Keep information per DPU
        dpu_info[k].rows_per_dpu = rows_per_dpu;
        dpu_info[k].rows_per_dpu_pad = rows_per_dpu_pad;
        dpu_info[k].prev_rows_dpu = prev_rows_dpu;
        dpu_info[k].nnz = nnz;

        // Find input arguments per DPU
        input_args[k].nrows = rows_per_dpu;
        input_args[k].tcols = s->ncols; 

#if BLNC_TSKLT_ROW
        // Load-balance rows across tasklets 
        partition_tsklt_by_row(part_info, rows_per_dpu, NR_TASKELTS);
#else
        // Load-balance nnzs across tasklets 
        partition_tsklt_by_nnz(s, part_info, rows_per_dpu, nnz, prev_rows_dpu, NR_TASKELTS);
#endif
        uint32_t t;
        for (t = 0; t < NR_TASKELTS; t++) {
            // Find input arguments per DPU
            input_args[k].start_row[t] = part_info->row_split_tasklet[t]; 
            input_args[k].rows_per_tasklet[t] = part_info->row_split_tasklet[t+1] - part_info->row_split_tasklet[t];

            if (input_args[k].rows_per_tasklet[t] > max_rows_per_tasklet)
                max_rows_per_tasklet = input_args[k].rows_per_tasklet[t];
        }
    }

    // Initializations for parallel transfers with padding needed
    if (max_rows_per_dpu % 2 != 0)
        max_rows_per_dpu++;
    if (max_nnz_ind_per_dpu % 2 != 0)
        max_nnz_ind_per_dpu++;
    if (max_nnz_val_per_dpu % (8 / byte_dt) != 0)
        max_nnz_val_per_dpu += ((8 / byte_dt) - (max_nnz_val_per_dpu % (8 / byte_dt)));
    if (max_rows_per_tasklet % (8 / byte_dt) != 0)
        max_rows_per_tasklet += ((8 / byte_dt) - (max_rows_per_tasklet % (8 / byte_dt)));

    // Re-allocations for padding needed
    s->rowptr = (uint32_t *) realloc(s->rowptr, (max_rows_per_dpu * nr_of_dpus * sizeof(uint32_t)));
    s->colind = (uint32_t *) realloc(s->colind, (max_nnz_ind_per_dpu * nr_of_dpus * sizeof(uint32_t)));
    s->values = (val_dt *) realloc(s->values, (max_nnz_val_per_dpu * nr_of_dpus * sizeof(val_dt)));
    y = (double*) malloc((uint64_t) ((uint64_t) nr_of_dpus * (uint64_t) max_rows_per_dpu) * (uint64_t) sizeof(val_dt)); 

    // Count total number of bytes to be transfered in MRAM of DPU
    unsigned long int total_bytes;
    total_bytes = ((max_rows_per_dpu) * sizeof(uint32_t)) + (max_nnz_ind_per_dpu * sizeof(uint32_t)) + (max_nnz_val_per_dpu * sizeof(val_dt)) + (s->ncols * sizeof(val_dt)) + (max_rows_per_dpu * sizeof(val_dt));
    assert(total_bytes <= DPU_CAPACITY && "Bytes needed exceeded MRAM size");

    // Copy input arguments to DPUs
	    k = 0;
    DPU_FOREACH(dpu_set, dpu, k) {
        input_args[k].max_rows = max_rows_per_dpu; 
        input_args[k].max_nnz_ind = max_nnz_ind_per_dpu; 
        DPU_ASSERT(dpu_prepare_xfer(dpu, input_args + k));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, "DPU_INPUT_ARGUMENTS", 0, sizeof(dpu_arguments_t), DPU_XFER_DEFAULT));


    // Copy input matrix to DPUs
	startTimer(&timer, 0);
    
	// Copy Rowptr 
	    k = 0;
    DPU_FOREACH(dpu_set, dpu, k) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, s->rowptr + dpu_info[k].prev_rows_dpu));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, (max_rows_per_dpu * sizeof(val_dt) + s->ncols * sizeof(val_dt)), max_rows_per_dpu * sizeof(uint32_t), DPU_XFER_DEFAULT));

    // Copy Colind
    k = 0;
    DPU_FOREACH(dpu_set, dpu, k) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, s->colind + s->rowptr[dpu_info[k].prev_rows_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, max_rows_per_dpu * sizeof(val_dt) + s->ncols * sizeof(val_dt) + max_rows_per_dpu * sizeof(uint32_t), max_nnz_ind_per_dpu * sizeof(uint32_t), DPU_XFER_DEFAULT));

    // Copy Values
    k = 0;
    DPU_FOREACH(dpu_set, dpu, k) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, s->values + s->rowptr[dpu_info[k].prev_rows_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, max_rows_per_dpu * sizeof(val_dt) + s->ncols * sizeof(val_dt) + max_rows_per_dpu * sizeof(uint32_t) + max_nnz_ind_per_dpu * sizeof(uint32_t), max_nnz_val_per_dpu * sizeof(val_dt), DPU_XFER_DEFAULT));
    stopTimer(&timer, 0);
	
	/*
	 * Print initial values of potential to output file.
	 */
	fprintf(output1, "0\t");
	for (i = 0; i < n; i++) {
		fprintf(output1, "%19.15f", u[i]);
	}	
		fprintf(output1, "\n");

	/*
	 * Temporal iteration.
	 */
	gettimeofday(&global_start, NULL);
	
	for (it = 0; it < total_time_steps; it++) {

	startTimer(&timer, 1);
    
	k = 0;
	//Send vector u to dpus
	DPU_FOREACH(dpu_set, dpu, k) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, u));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, max_rows_per_dpu * sizeof(val_dt), s->ncols * sizeof(val_dt), DPU_XFER_DEFAULT));
    stopTimer(&timer, 1);
	
	average_vector_allocation_time += timer.time[1];
	//printf("Vector Allocation in each iteration ");
	//printTimer(&timer, 1);

	//Run kernel on DPUs
   startTimer(&timer, 2);
    DPU_ASSERT(dpu_launch(dpu_set, DPU_SYNCHRONOUS));
    stopTimer(&timer, 2);
	
	average_run_time += timer.time[2];
	//printf("Kernel in a signle iteration ");
	//printTimer(&timer, 2);
	
	//Gather y output vector
    startTimer(&timer, 3);
    k = 0;
    DPU_FOREACH(dpu_set, dpu, k) {
        DPU_ASSERT(dpu_prepare_xfer(dpu, y + (k * max_rows_per_dpu)));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_FROM_DPU, DPU_MRAM_HEAP_POINTER_NAME, 0, max_rows_per_dpu * sizeof(val_dt), DPU_XFER_DEFAULT));
    stopTimer(&timer, 3);
	average_output_gather_time += timer.time[3];
	
	//printf("Gather output in a signle iteration ");
	//printTimer(&timer, 3);
	//printf("\n");


		i=0;
		 for(int n_dp=0; n_dp<nr_of_dpus; n_dp++){
			cur_rows = dpu_info[n_dp].rows_per_dpu;
			for(int z=0; z<cur_rows; z++){
			u_next[i] = u[i] + dt * (mu - u[i]) + dt * (y[n_dp * max_rows_per_dpu + z] - sum_s[i] * u[i]) / divide[i];

			/*
			 * Update network elements and set u[i] = 0 if u[i] > uth
			 */
			if (u_next[i] > uth) {
				u_next[i] = 0.0;
				/*
				 * Calculate omega's.
				 */
				if (it >= ttransient) {
					omega1[i] += 1.0;
				}
			}
			i++;
			}
		}

		 // Exchange u and u_next
		 
		temp = u;
		u = u_next;
		u_next = temp;

		/*
		 * Print out of results.
		 */
		
		if ((it + 1) % sampling == 0) {
			
			total_average_run_time += average_run_time;
			total_vector_allocation_time += average_vector_allocation_time;
			total_output_gather_time += average_output_gather_time;
			
			//printf("\nTime is %ld: \n", it + 1);
			/*
			
			THESE PRINTF ARE IN CASE SOMEONE WHANTS TO PRINT THE AVERAGE TIME IN EVERY 100 TIME STEPS (FOR EVERY SAMPLE)

			printf("The average vector allocation time is (msec): %lf\n", (average_vector_allocation_time / sampling)*0.001);
			printf("The average kernel run time is (sec): %lf\n", (average_run_time / sampling)*0.000001);
			printf("The output gather run time is (msec): %lf\n\n", (average_output_gather_time / sampling)*0.001);
			
			*/

			average_run_time = 0;
			average_vector_allocation_time = 0;
			average_output_gather_time = 0;

			gettimeofday(&IO_start, NULL);
			fprintf(output1, "%ld\t", it + 1);
			for (i = 0; i < n; i++) {
				fprintf(output1, "%19.15f", u[i]);
			}
			fprintf(output1, "\n");						

			if (it > ttransient) {
				time = (double)it * dt;
				fprintf(output2, "%ld\t", it + 1);
				for (i = 0; i < n; i++) {
					omega[i] = 2.0 * M_PI * omega1[i] / (time - ttransient * dt);
					fprintf(output2, "%19.15f", omega[i]);
				}
				fprintf(output2, "\n");
			}
			gettimeofday(&IO_end, NULL);
			IO_usec += ((IO_end.tv_sec - IO_start.tv_sec) * 1000000.0 + (IO_end.tv_usec - IO_start.tv_usec));
		}
	}
	
	printf("Total average vector allocation time (msec) %lf\n", (total_vector_allocation_time/ total_time_steps)*0.001);
	printf("Total average Kernel run time (sec) %lf\n", (total_average_run_time / total_time_steps)*0.000001);
	printf("Total average output run time (msec) %lf\n", (total_output_gather_time/total_time_steps)*0.001);

	printf("\nTotal time for alloactions: (sec): %lf", total_vector_allocation_time*0.000001);
	printf("\nTotal time for Kernel: (sec): %lf ", total_average_run_time*0.000001);
	printf("\nTotal time for Gather: (sec): %lf \n", total_output_gather_time*0.000001);
	
	printf("\nMatrix allcation to dpu (msec): ");
	printf("%lf",timer.time[0]*0.001);
	printf("\n\n");

	gettimeofday(&global_end, NULL);
	global_usec = ((global_end.tv_sec - global_start.tv_sec) * 1000000.0 + (global_end.tv_usec - global_start.tv_usec));

	printf("Time for calculations = %13.6f sec\n", (global_usec - IO_usec) / 1000000.0);
	printf("Time for I/O          = %13.6f sec\n", IO_usec / 1000000.0);
	printf("Total execution time  = %13.6f sec\n", global_usec / 1000000.0);

	fclose(output1);
	fclose(output2);

	return 0;
}