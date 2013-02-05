#include <mpi.h>
#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <stdlib.h>

#define NPROCS 6
#define LONG_SIDE 128	// size of coordinates
#define SHORT_SIDE 64
#define FIELD 0	// field process has rank of 0
#define INVALID -1
#define NO_OF_ROUNDS 900
#define MAX_STEPS 10


typedef struct {

	int ball_coords[2] ;
	int dist ;
	int final_coords[2] ;
	int initial_coords[2] ;
	int passed ;
	int rank ;
	int reached ;
	int win ;

} Player;

void init_coords(int coords[]) ;
void init_player(Player* process) ;
int player_run(Player* runner) ;
void print_stats(Player player_info[], int win) ;
MPI_Datatype createMPIStruct() ;
void inform_field(Player* process) ;

int main(int argc, char *argv[]) {

	int i, numtasks, win, rounds ;
	Player process ;
	Player player_info[NPROCS] ;
	MPI_Datatype player_MPI ;


	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	if (numtasks!=NPROCS) {
		printf("Must specify %d processors. Terminating.\n",NPROCS);
		MPI_Finalize();
		return 1;	
	}

	player_MPI = createMPIStruct() ;
	MPI_Type_commit(&player_MPI); 

	init_player(&process) ;

	for (rounds=1; rounds<=NO_OF_ROUNDS; rounds++){

		// Broadcast ball location
		MPI_Bcast(process.ball_coords, 2, MPI_INT, FIELD, MPI_COMM_WORLD); 
		player_run(&process) ; 

		// Players report their coordinates
		MPI_Gather(&process, 1, player_MPI, &player_info, 1, player_MPI, FIELD, MPI_COMM_WORLD); 
		process.win = decide_winner(player_info); 

		// Field tells everyone who is the winner
		MPI_Bcast(&process.win, 2, MPI_INT, FIELD, MPI_COMM_WORLD);

		// Winner passes the ball
		if (process.rank == process.win) {
			init_coords(process.ball_coords);
			process.passed++ ;			
		}

		inform_field(&process);

		if (process.rank == FIELD) {
			printf("%d\n", rounds);
			printf("%d %d\n", process.ball_coords[0], process.ball_coords[1]);
			print_stats(player_info, process.win);
		}

		// new coords become the old coords of next round
		else {
			process.initial_coords[0] = process.final_coords[0] ;
			process.initial_coords[1] = process.final_coords[1] ;

		}

	}

	MPI_Finalize();
	return 0;
}


// initializes player information
void init_player(Player* process) {
	int rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	process->rank = rank;

	// use rank to make the seed random for each process
	srand(time(NULL) * rank);

	if (rank == FIELD) {
		process->initial_coords[0] = INVALID ;
		process->initial_coords[1] = INVALID ;
		process->final_coords[0] = INVALID ;
		process->final_coords[1] = INVALID ;
		process->dist = INVALID;
		init_coords(process->ball_coords) ;	
		process->reached = INVALID ;
		process->passed = INVALID ;
		process->win = INVALID ;
	}
	else
	{
		init_coords(process->initial_coords) ;
		process->final_coords[0] = INVALID ;
		process->final_coords[1] = INVALID ;
		process->dist = 0 ;
		process->ball_coords[0] = INVALID ;
		process->ball_coords[1] = INVALID ;
		process->reached = 0 ;
		process->passed = 0 ;
		process->win = INVALID ;
	}

}

void init_coords(int coords[])
{
	coords[0] = rand() % (LONG_SIDE + 1) ;
	coords[1] = rand() % (SHORT_SIDE + 1) ;
}

// Updates coordinates of player after running
int player_run(Player* runner)
{	
	int steps_to_run = (rand() % MAX_STEPS) + 1;
	int x_diff = runner->ball_coords[0] - runner->initial_coords[0] ;
	int y_diff = runner->ball_coords[1] - runner->initial_coords[1] ;
	int steps_to_ball = abs(x_diff) + abs(y_diff) ;

	// if can reach ball, run to ball position
	if (steps_to_ball <= steps_to_run){
		runner->final_coords[0] = runner->ball_coords[0] ;
		runner->final_coords[1] = runner->ball_coords[1] ;
		runner->reached++ ;
		steps_to_run = steps_to_ball;
	}

	else	// run the shortest path towards ball
	{
		int x_run = (float)x_diff/steps_to_ball * steps_to_run ;
		int y_run = steps_to_run - abs(x_run) ;

		// run in right direction
		if (y_diff<0)
			y_run *=-1 ;

		runner->final_coords[0] = runner->initial_coords[0] + x_run ;
		runner->final_coords[1] = runner->initial_coords[1] + y_run ;
	}

	runner->dist += steps_to_run ;
	return steps_to_run;
}

// Decides who gets possession of the ball
int decide_winner(Player player_info[])
{
	int players_at_ball[NPROCS] = {0}, count, i ;

	for (i=0, count=0; i < NPROCS; i++)
	{

		if (player_info[i].final_coords[0] == player_info[i].ball_coords[0] &&
			player_info[i].final_coords[1] == player_info[i].ball_coords[1] &&
			player_info[i].rank!=FIELD){
				players_at_ball[count] = i ;
				count++;
		}

	}

	if (count > 0) {
		int winner_index = players_at_ball[rand() % count];
		player_info[winner_index].passed ++ ;
		return player_info[winner_index].rank ;
	}
	return -1 ;
}

// print player stats every round
void print_stats(Player player_info[], int win)
{
	int i, at_ball, pass_ball ;
	for (i=0; i < NPROCS; i++)
	{
		if (player_info[i].rank != FIELD) {
			if (player_info[i].final_coords[0] == player_info[i].ball_coords[0] &&
				player_info[i].final_coords[1] == player_info[i].ball_coords[1]) {
					at_ball = 1;
					if (player_info[i].rank==win)
						pass_ball = 1;
					else
						pass_ball = 0;
			}
			else {
				at_ball = 0;
				pass_ball = 0;
			}


			printf("%d %d %d ", player_info[i].rank-1, player_info[i].initial_coords[0], player_info[i].initial_coords[1]);
			printf("%d %d %d %d ", player_info[i].final_coords[0], player_info[i].final_coords[1], at_ball, pass_ball);
			printf("%d %d %d\n", player_info[i].dist, player_info[i].reached, player_info[i].passed); 
		}
	}

}

// define a MPI struct for communication
MPI_Datatype createMPIStruct(){
	int blocklen[1] = {11} ;
	MPI_Datatype oldtype[1] = {MPI_INT}, newtype ;
	MPI_Aint disp[1]; 	
	disp[0] = 0;

	MPI_Type_struct(1, blocklen, disp, oldtype, &newtype); 

	return newtype;

}

// winner scatters new ball position to all processes
void inform_field(Player* process)
{
	int ball_info[NPROCS][2], i ;
	if (process->rank == process->win)
	{
		for (i=0; i<NPROCS; i++) {
			ball_info[i][0] = -1;
			ball_info[i][1] = -1;
		}

		ball_info[FIELD][0] = process->ball_coords[0];
		ball_info[FIELD][1] = process->ball_coords[1];
		ball_info[process->win][0] = process->ball_coords[0];
		ball_info[process->win][1] = process->ball_coords[1];

	}
	if (process->win!=-1)	// only scatter if there is winner
		MPI_Scatter (ball_info, 2, MPI_INT, process->ball_coords, 2, MPI_INT, process->win, MPI_COMM_WORLD) ;
}


