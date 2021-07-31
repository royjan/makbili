#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "myProto.h"
#include <omp.h>
#include <stddef.h>

#include <unistd.h>

//my
int new_main(int argc, char *argv[]);
void new_slave(char *seq1, char *seq2);
void new_master(char **seq1, char **seq2, char **procedure);
void checkOneCase(int my_rank, int batch, int x, Scores *scores_each_t,
		char *seq1, char *seq2, int size_seq2);
//

void createScoreType(MPI_Datatype *scoreMpiType);
void createWeightType(MPI_Datatype *weightMpiType);

void readFromFile(const char *file_path, Weights *weights, char **seq1,
		char **seq2, char **procedure, int *size_seq1, int *size_seq2) {
	int i;
	int size_proc;
	FILE *file = fopen(file_path, "r");
	if (!file) {
	}
	i = fscanf(file, "%lf %lf %lf %lf", &weights->w1, &weights->w2,
			&weights->w3, &weights->w4);
	*seq1 = (char*) malloc(sizeof(char) * MAX_CHARS_FOR_SEQ1);
	i = fscanf(file, "%s", *seq1);
	*size_seq1 = strlen(*seq1);

	if (*size_seq1 > MAX_CHARS_FOR_SEQ1) {
		printf(
				"\n****ERROR! cannot read the file, the string for Sequence 1 is too long of length:%d!****\n",
				*size_seq1);
		exit(1);
	}
	*seq2 = (char*) malloc(sizeof(char) * MAX_CHARS_FOR_SEQ2);
	i = fscanf(file, "%s", *seq2);
	*size_seq2 = strlen(*seq2);
	if (*size_seq1 > MAX_CHARS_FOR_SEQ2) {
		printf(
				"\n****ERROR! cannot read the file, the string for Sequence 2 is too long of length:%d!****\n",
				*size_seq2);
		exit(1);
	}

	*procedure = (char*) malloc(sizeof(char) * MAX_CHARS_FOR_PROCEDURE);
	i = fscanf(file, "%s", *procedure);
	size_proc = strlen(*procedure);
	if (size_proc > MAX_CHARS_FOR_PROCEDURE) {
		printf(
				"\n****ERROR! cannot read the file, the string for procedure is too long of length:%d!****\n",
				size_proc);
		exit(1);
	}
	fclose(file);
}

void writeToFile(FILE *file, Scores *scores, char *seq2) {
	fprintf(file,
			"The mutant is : %s\t|Best Score: %1.3f\t|Best offset: %d\t\n",
			seq2, scores->score, scores->offset);

}

int checkSemiOrSemiConservativeGroups1(char c1, char c2, const char *group[],
		int size) {
	int i, j;
	char found_c1;
	char found_c2;
	int size_g;

	for (i = 0; i < size; i++) {
		found_c1 = '-';
		found_c2 = '-';
		size_g = strlen((char*) group[i]);

		for (j = 0; j < size_g; j++) {
			if (found_c1 == '-' && c1 == group[i][j]) {
				found_c1 = c1;
			}
			if (found_c2 == '-' && c2 == group[i][j]) {
				found_c2 = c2;
			}
			if (found_c1 != '-' && found_c2 != '-') {
				return 1;

			}
		}
	}
	return 0;
}

int calcMaxOffsetForEachSeq(char *seq1, char *seq2) {
	int size_seq1;
	int size_seq2;
	int maxOffsetEachSeq2;
	size_seq1 = strlen(seq1);
	size_seq2 = strlen(seq2);
	maxOffsetEachSeq2 = abs(size_seq1 - size_seq2);
	return maxOffsetEachSeq2;
}

void addMutant(char *seq1, char *seq2, int size_seq2, char *new_seq2) {

	const char *cons_Groupss[9] = { "NDEQ", "MILV", "FYM", "NEQK", "QHRK", "HY",
			"STA", "NHQK", "MILF" };
	int j, k, i, d;
	int seq2_flag = 0, seq1_flag = 0, flag = 1;

	//printf("%lu\n",sizeof(cons_Groupss)/sizeof(char*));
	for (j = 0; j < size_seq2; j++) {

		flag = 1;
		if (seq1[j] != seq2[j]) {

			for (d = 0; d < sizeof(cons_Groupss) / sizeof(char*); d++) {
				seq2_flag = 0, seq1_flag = 0;

				for (i = 0; i < strlen(cons_Groupss[d]); i++) {

					if (seq2[j] == cons_Groupss[d][i]) {
						seq2_flag = 1;

						if (seq1_flag == 1) {
							seq1_flag = 0;
							seq2_flag = 0;
							flag = 0;
							printf("they are both under the same group\n");

						}

					}
					if (seq1[j] == cons_Groupss[d][i]) {
						seq1_flag = 1;

						if (seq2_flag == 1) {
							seq1_flag = 0;
							seq2_flag = 0;
							flag = 0;
							printf("they are both under the same group\n");

						}

					}
				}

			}
			if ((!(seq2_flag == 1 && seq1_flag == 1)) && flag == 1) {
				//printf("%c",seq2[j]);
				//printf("%c",seq1[j]);
				seq2[j] = seq1[j];

			}

		}
	}
	//printf("%s\n",seq2);
	for (j = 0; j < size_seq2 + 1; j++) {
		new_seq2[j] = seq2[j];
	}

}

void addMutantsEachSeq(char *seq2, int size_seq2, char **new_seq2, int *i,
		Scores **scores) {

	int j, k;
	int size_new_seq = size_seq2 + 1;
	*new_seq2 = (char*) calloc(size_new_seq, sizeof(char));
	for (j = 0, k = 0; j < size_seq2 + 1; j++, k++) {
		if (*i == j) {
			(*new_seq2)[*i] = '-';
			(*scores)->mute_loc = *i;
			j += 1;
		}
		(*new_seq2)[j] = seq2[k];
	}
}

double compareDiffrentLenSequences(Weights weights, char *seq1, char *seq2,
		int size_seq2, Scores **scores) {
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	int maxOffset;
	int i;
	int tid = 0;
	char *sub;
	double max_score = -INFINITY;

	double *max_scores_each_t;
	int *offset_each_t;
	int n_ot = 4;
	omp_set_num_threads(n_ot);
	max_scores_each_t = (double*) calloc(n_ot, sizeof(double));
	offset_each_t = (int*) calloc(n_ot, sizeof(int));
	(*scores)->offset = 0;
	maxOffset = calcMaxOffsetForEachSeq(seq1, seq2);

#pragma omp parallel private(i, tid, sub)
	{
		sub = (char*) malloc(size_seq2 * sizeof(char));
		tid = omp_get_thread_num();
		max_scores_each_t[tid] = -INFINITY;

#pragma omp for
		for (i = 0; i < maxOffset + 1; i++) {
			double score = 0.0;
			substring(seq1, sub, i + 1, size_seq2);
			score = 2.0;//compareSameLenSequences(weights, sub, seq2, size_seq2);
			if (score > max_scores_each_t[tid]) {
				max_scores_each_t[tid] = score;
				offset_each_t[tid] = i;
			}
		}
	}

	for (i = 0; i < n_ot; i++) {
		if (max_scores_each_t[i] > max_score) {
			max_score = max_scores_each_t[i];
			(*scores)->offset = offset_each_t[i];
		}
	}

	return max_score;
}

double compareSameLenSequencesNew(char *seq1, char *seq2, int size_seq2,
		int offset) {
	const char *cons_Groups[9] = { "NDEQ", "MILV", "FYM", "NEQK", "QHRK", "HY",
			"STA", "NHQK", "MILF" };
	const char *semi_Cons_Groups[11] = { "SAG", "SGND", "NEQHRK", "ATV", "STPA",
			"NDEQHK", "HFY", "CSA", "STNK", "SNDEQK", "FVLIM" };
	int i;
	Counts counts;
	double score = 0.0;
	counts.countStars = 0;
	counts.countDots = 0;
	counts.countColons = 0;
	counts.countSpaces = 0;
	Weights weights;
	weights.w1 = 1;
	weights.w2 = 1;
	weights.w3 = 1;
	weights.w4 = 1;
	int current_offset = offset;
	printf("\n\n");
	for (i = 0; i < size_seq2; i++) {
		current_offset = i + offset;
		printf("Long:%c Short:%c\n", seq1[current_offset], seq2[i]);
		if (seq1[current_offset] == seq2[i]) {
			counts.countStars++;
		} else if (checkSemiOrSemiConservativeGroups1(seq1[current_offset],
				seq2[i], cons_Groups, 9) == 1) {
			counts.countColons++;
		} else if (checkSemiOrSemiConservativeGroups1(seq1[current_offset],
				seq2[i], semi_Cons_Groups, 11) == 1) {
			counts.countDots++;
		} else {
			counts.countSpaces++;
		}
	}
	printf("w1=%1.3f | w2=%1.3f | w3=%1.3f | w4=%1.3f\n", weights.w1,
			weights.w2, weights.w3, weights.w4);
	printf("starts:%d | colons:%d | dots:%d | spaces:%d\n", counts.countStars,
			counts.countColons, counts.countDots, counts.countSpaces);
	score = (weights.w1 * counts.countStars) - (weights.w2 * counts.countColons)
			- (weights.w3 * counts.countDots)
			- (weights.w4 * counts.countSpaces);
	return score;
}

void substring(char *s, char *sub, int p, int l) {
	int c = 0;

	while (c < l) {
		sub[c] = s[p + c - 1];
		c++;
	}
	sub[c] = '\0';
}

void createWeightType(MPI_Datatype *weightMpiType) {
	int blocklengths[WEIGHT_NUM_ATTRIBUTES] = WEIGHT_BLOCK_LENGTH;
	MPI_Datatype types[WEIGHT_NUM_ATTRIBUTES] = WEIGHT_TYPE;
	MPI_Aint offsets[WEIGHT_NUM_ATTRIBUTES];
	offsets[0] = offsetof(Weights, w1);
	offsets[1] = offsetof(Weights, w2);
	offsets[2] = offsetof(Weights, w3);
	offsets[3] = offsetof(Weights, w4);
	MPI_Type_create_struct(WEIGHT_NUM_ATTRIBUTES, blocklengths, offsets, types,
			weightMpiType);
	MPI_Type_commit(weightMpiType);
}

void createScoreType(MPI_Datatype *scoreMpiType) {
	int blocklengths[SCORE_NUM_ATTRIBUTES] = SCORE_BLOCK_LENGTH;
	MPI_Datatype types[SCORE_NUM_ATTRIBUTES] = SCORE_TYPE;
	MPI_Aint offsets[SCORE_NUM_ATTRIBUTES];
	offsets[0] = offsetof(Scores, score);
	offsets[1] = offsetof(Scores, mute_loc);
	offsets[2] = offsetof(Scores, offset);
	offsets[3] = offsetof(Scores, num_of_Seq);
	MPI_Type_create_struct(SCORE_NUM_ATTRIBUTES, blocklengths, offsets, types,
			scoreMpiType);
	MPI_Type_commit(scoreMpiType);
}

int main(int argc, char *argv[]) {

	int my_rank;
	int p;
	char *seq1;
	//int numOfSeq2;
	char *procedure;
	char *seq2;
	char *one_seq = NULL;
	double t1, t2, time;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Comm_size(MPI_COMM_WORLD, &p);

	t1 = MPI_Wtime();

	if (my_rank != 0) {

		new_slave(seq1, seq2);

	} else {

		new_master(&seq1, &seq2, &procedure);
		t2 = MPI_Wtime();
		time = (t2 - t1) / 60; // get the execution time in minutes
		printf("Executation time = %1.3f m \n", time);
	}

	MPI_Finalize();

	return 0;
}

void new_master(char **seq1, char **seq2, char **procedure) {
	int i;
	int size;
	MPI_Datatype weightMPIType;
	createWeightType(&weightMPIType);
	int size_seq2;
	int size_seq1;
	int batch = 0;

	MPI_Status status;
	MPI_Datatype scoreMPIType;
	createScoreType(&scoreMPIType);

	Scores scores;
	Weights weights;
	Counts counts;

	readFromFile("input.txt", &weights, seq1, seq2, procedure, &size_seq1,
			&size_seq2);
	int p;
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	int maxOffset = calcMaxOffsetForEachSeq(*seq1, *seq2);
	batch = ceil(maxOffset / (p - 1));

	int task_count = 0;
	int term_count = 0;
	int TERMINATION_TAG = 1;

	MPI_Bcast(&size_seq1, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Bcast(*seq1, size_seq1, MPI_CHAR, 0, MPI_COMM_WORLD);

	MPI_Bcast(&size_seq2, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(*seq2, size_seq2, MPI_CHAR, 0, MPI_COMM_WORLD);

	MPI_Bcast(&weights, 1, weightMPIType, 0, MPI_COMM_WORLD);

	MPI_Bcast(&batch, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//output file 
	const char *file_path;
	file_path = "output.txt";
	FILE *file = fopen(file_path, "w");

	//int num_seq;
	writeToFile(file, &scores, *seq1);

	Scores *scores_arr;
	//scores_arr = (Scores*) malloc(*numOfSeq2 * sizeof(Scores));

	fclose(file);

}

void new_slave(char *seq1, char *seq2) {

	int dest, i;
	dest = 0;
	int f;
	int my_rank;
	int size_seq1;
	int size_seq2;
	int batch = 0;

	Weights weights;
	Counts counts;

	MPI_Datatype scoreMPIType;
	createScoreType(&scoreMPIType);
	MPI_Datatype weightMPIType;
	createWeightType(&weightMPIType);

	int termination_tag = 1;

	MPI_Status status;

	int tag = 0;
	int size;

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Bcast(&size_seq1, 1, MPI_INT, 0, MPI_COMM_WORLD);
	seq1 = (char*) calloc(size_seq1, sizeof(char));
	MPI_Bcast(seq1, size_seq1, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(&size_seq2, 1, MPI_INT, 0, MPI_COMM_WORLD);
	seq2 = (char*) calloc(size_seq2, sizeof(char));
	MPI_Bcast(seq2, size_seq2, MPI_CHAR, 0, MPI_COMM_WORLD);

	MPI_Bcast(&weights, 1, weightMPIType, 0, MPI_COMM_WORLD);
	MPI_Bcast(&batch, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//from here its thread shit 
	int tid = 0;
	char *sub;

	int n_ot = 1;
	omp_set_num_threads(n_ot);

	Scores *scores_each_t;
	scores_each_t = (Scores*) calloc(batch, sizeof(Scores));
	int offset;

#pragma omp parallel
	{

#pragma omp single
		{
			for (offset = (my_rank - 1) * batch; offset <= my_rank * batch;
					offset++) {
				int d = offset; //dont touch!
#pragma omp task
				{
					checkOneCase(my_rank, batch, d, scores_each_t, seq1, seq2,
							size_seq2);
				}
			}

		}
	}

	Scores tempMax;
	tempMax.score = -INFINITY;
	tempMax.offset = -1;
	//here find max in result arr
	for (f = 0; f < batch; f++) {
		if (scores_each_t[f].score > tempMax.score) {
			tempMax = scores_each_t[f];
		}
	}

	printf("MPI%2d|Max offset %3d|with score %1.3f \n", my_rank, tempMax.offset,
			tempMax.score);
	exit(0);

}

void checkOneCase(int my_rank, int batch, int x, Scores *scores_each_t,
		char *seq1, char *seq2, int size_seq2) {
	double score = compareSameLenSequencesNew(seq1, seq2, size_seq2, x);
	printf("SCORE: %1.3f\n", score);
	scores_each_t[x % batch].score = score;
	scores_each_t[x % batch].offset = x;
}

