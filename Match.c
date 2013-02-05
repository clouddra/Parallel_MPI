#include <mpi.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
//#include <sys/time.h> //Timing Code

#define NPROCS 12
#define LONG_SIDE 128	// size of coordinates
#define SHORT_SIDE 64
#define RANK 1	// size of rank
#define FIELD_0 0	// field process has rank of 0
#define FIELD_1 1	// field process has rank of 1
#define FIELD 0		
#define INVALID -1	
#define INVALID_COORDS -1000	// Invalid value for coordinates
#define PASS_FACTOR 1	
#define HALF_MATCH 2700

#define T_inform_ball_coords 1
#define T_report_pos 2
#define T_inform_winner 3
#define T_players_on_field 4
#define T_inform_shoot_pos 5
#define T_update_field0_player 6
#define T_winner_field1 7
#define T_update_field0_ball 8

//Timing Code
/*
long long wall_clock_time()
{
#ifdef __linux__
struct timespec tp;
clock_gettime(CLOCK_REALTIME, &tp);
return (long long)(tp.tv_nsec + (long long)tp.tv_sec * 1000000000ll);
#else
#warning "Your timer resoultion might be too low. Compile on Linux and link with librt"
struct timeval tv;
gettimeofday(&tv, NULL);
return (long long)(tv.tv_usec * 1000 + (long long)tv.tv_sec * 1000000000ll);
#endif
}
*/

typedef struct {

	int global_rank;
	int player_id;
	int initial_coords[2] ;
	int final_coords[2] ;
	int reached ;
	int challenge ;
	int ball_coords[2] ;	// shoot location for players, ball location for field
	int team;

} Player_Match;

typedef struct {

	int global_rank;
	int court_rank;	
	int speed;
	int dribble;
	int shooting;
	int post[2];
	MPI_Comm team_comm; 
	MPI_Comm court_comm;
	int court_group_size;


} Player_Stats;


void init_coords(int coords[])
{
	coords[0] = rand() % (LONG_SIDE + 1) ;
	coords[1] = rand() % (SHORT_SIDE + 1) ;
}

// initializa player stats
void init_player_stats(Player_Stats* process) {

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// ensure that the seed of each rank is different
	srand(time(NULL) * rank+1);
	process->global_rank = rank;

	// field does not have stats
	if (rank == FIELD_0 || rank == FIELD_1) {
		process->dribble = INVALID;
		process->speed = INVALID ;
		process->shooting = INVALID ;


	}
	else
	{
		process->dribble = rand() % 5 + 1;	// 1-5
		process->speed = 5 + rand()%5 ;		// 5-9
		process->shooting = 15- process->dribble - process->speed;	// 1-9
	}

}

// initialize player's match information (coords)
void init_player_match(Player_Match* process) {
	int rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	process->global_rank = rank;

	// field does not have coordinates
	if (rank == FIELD_0 || rank == FIELD_1) {
		process->initial_coords[0] = INVALID_COORDS ;
		process->initial_coords[1] = INVALID_COORDS ;
		process->final_coords[0] = INVALID_COORDS ;
		process->final_coords[1] = INVALID_COORDS ;
		process->ball_coords[0] = LONG_SIDE/2 ;
		process->ball_coords[1] = SHORT_SIDE/2 ;
		process->reached = INVALID ;
		process->challenge = INVALID ;

	}
	else
	{
		// initial coordinates of players are random
		init_coords(process->initial_coords) ;
		process->final_coords[0] = INVALID_COORDS ;
		process->final_coords[1] = INVALID_COORDS ;
		process->ball_coords[0] = INVALID_COORDS ;
		process->ball_coords[1] = INVALID_COORDS ;
		process->reached = 0 ;
		process->challenge = 0 ;

	}

}

// initialize the information about the player's team
void init_team(Player_Stats* stats, Player_Match* match)
{

	if (stats->global_rank == FIELD_0 || stats->global_rank == FIELD_1)
		match->team = 0;
	else if ((stats->global_rank) <= NPROCS/2) {
		match->team = 1;
		stats->post[0] = 0;
		stats->post[1] = SHORT_SIDE/2;
	}
	else {
		match->team  = 2;
		stats->post[0] = 128;
		stats->post[1] = SHORT_SIDE/2;
	}
	MPI_Comm_split(MPI_COMM_WORLD, match->team, stats->global_rank, &(stats->team_comm));
	MPI_Comm_rank(stats->team_comm, &(match->player_id));


}

// sets up an wall around the court. Ensures that coodinates does not go beyond boundaries
void limit_coords(int coords[])
{
	coords[0] =  (LONG_SIDE < coords[0]) ? LONG_SIDE : coords[0];
	coords[0] =  (0 > coords[0]) ? 0 : coords[0];
	coords[1] =  (SHORT_SIDE < coords[1]) ? SHORT_SIDE : coords[1];
	coords[1] =  (0 > coords[1]) ? 0 : coords[1];

}


void setup_field_comm(Player_Stats* stats, Player_Match* match, MPI_Comm* comm)
{
	int field_number, long_side ;
	long_side = match->final_coords[0];


	if (long_side==INVALID_COORDS)
		long_side = match->initial_coords[0] ;

	if (stats->global_rank == FIELD_0 || stats->global_rank == FIELD_1)
		field_number = stats->global_rank;
	else if (long_side < LONG_SIDE/2)
		field_number = FIELD_0;
	else
		field_number = FIELD_1;


	MPI_Comm_split(MPI_COMM_WORLD, field_number, stats->global_rank, comm);
	MPI_Comm_rank(*comm, &(stats->court_rank));
	MPI_Comm_size(*comm, &(stats->court_group_size));

}

// Updates coordinates of player after running
int player_run(Player_Match* runner, Player_Stats* stats)
{	
	int steps_to_run = stats->speed;
	int x_diff = runner->ball_coords[0] - runner->initial_coords[0] ;
	int y_diff = runner->ball_coords[1] - runner->initial_coords[1] ;
	int steps_to_ball = abs(x_diff) + abs(y_diff) ;

	// if can reach ball, run to ball position
	if (steps_to_ball <= steps_to_run){
		runner->final_coords[0] = runner->ball_coords[0] ;
		runner->final_coords[1] = runner->ball_coords[1] ;
		runner->reached = 1 ;
		runner->challenge = (rand() % 10 + 1) *  stats->dribble ;
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
		runner->reached = 0 ;
		runner->challenge = -1 ;

	}

	runner->ball_coords[0] = INVALID_COORDS ;
	runner->ball_coords[1]= INVALID_COORDS;
	return steps_to_run;
}

// calculates the chance of ball landing in target position
float shoot_chance(int shooter_pos[], int target[], int skill){
	float prob;
	int d;
	d = abs(shooter_pos[0]-target[0]) + abs(shooter_pos[1]-target[1]) ;
	//printf("dist %d\n", d);

	if (d==0)
		return 100.0;
	else {
		prob = (10 + (float)90*skill)/(0.5*d*sqrt((float)d)-0.5);
		prob =  (100.0 < prob) ? 100.0 : prob;
		//printf("prob %f\n", prob);
		return prob;
	}
}

// define a MPI struct for communication
MPI_Datatype createPlayerMatchStruct(){
	int blocklen[1] = {11} ;
	MPI_Datatype oldtype[1] = {MPI_INT}, newtype ;
	MPI_Aint disp[1]; 	
	disp[0] = 0;

	MPI_Type_struct(1, blocklen, disp, oldtype, &newtype); 

	return newtype;

}

// Decides who gets possession of the ball
int decide_winner(Player_Match match_info[], int group_size)
{
	int max_challenge, winner_index, count, i ;
	int players_at_ball[NPROCS-2] ;

	for (i=0, count=0, max_challenge=-1; i < group_size; i++)
	{
		if (match_info[i].challenge > max_challenge) {
			players_at_ball[0] = i;

			max_challenge = match_info[i].challenge;
			count=1;


		}

		else if (match_info[i].challenge == max_challenge)
		{
			players_at_ball[count] = i;
			count++;
		}
	}

	// player rank starts from 1 since field process is not in match_info
	if (max_challenge > 0) {
		winner_index = rand()%count ;

		return players_at_ball[winner_index] + 1;	
	}

	return INVALID;
}

// Prints information about each round
int print_round_info(Player_Match match_info[], int score[], int final_ball_coords[])
{
	int i, sort_index, scored;
	Player_Match sorted_info[NPROCS-2];
	int team_offset, shoot_ball;

	// sort the players by team and id
	for (i=0, scored=0; i<NPROCS-2; i++) {
		sort_index = (match_info[i].team - 1) * (NPROCS-2)/2 + match_info[i].player_id;
		sorted_info[sort_index] = match_info[i];

		if (match_info[i].ball_coords[1]!=INVALID_COORDS) {
			// the ball coords of the player is one of the posts
			if ((match_info[i].ball_coords[0] == 0 || match_info[i].ball_coords[0] == LONG_SIDE) && match_info[i].ball_coords[1] == SHORT_SIDE/2) {
				scored++;
				// update score based on distance
				if (abs(match_info[i].ball_coords[0] - match_info[i].final_coords[0]) + abs(match_info[i].ball_coords[1] - match_info[i].final_coords[1]) > 24)
					score[match_info[i].team]+=3 ;
				else 
					score[match_info[i].team]+=2 ;
			}
			else {
				sorted_info[sort_index].ball_coords[0] = abs(match_info[i].ball_coords[0]) ;
				sorted_info[sort_index].ball_coords[1] = abs(match_info[i].ball_coords[1]) ;
			}
		}

		else 
		{
			sorted_info[sort_index].ball_coords[0] = -1 ;
			sorted_info[sort_index].ball_coords[1] = -1 ;
		}
	}

	printf("%d %d\n", score[1], score[2]) ;
	printf("%d %d\n", final_ball_coords[0], final_ball_coords[1]) ;

	for (i=0; i<NPROCS-2; i++) {

		printf("%d ", sorted_info[i].player_id); 
		printf("%d %d ", sorted_info[i].initial_coords[0], sorted_info[i].initial_coords[1]);
		printf("%d %d ", sorted_info[i].final_coords[0], sorted_info[i].final_coords[1]);
		printf("%d ", sorted_info[i].reached);
		if (sorted_info[i].ball_coords[0] == -1)
			printf("0 ");
		else 
			printf("1 ");
		printf("%d ", sorted_info[i].challenge);
		printf("%d %d\n", sorted_info[i].ball_coords[0], sorted_info[i].ball_coords[1]);

	}

	// reset ball position if someone scored
	if (scored!=0) {
		final_ball_coords[0] = LONG_SIDE/2 ;
		final_ball_coords[1] = SHORT_SIDE/2 ;
	}
}

// determine target position and whether the ball landed in target position
void shoot(Player_Match other_players[], Player_Match* shooter, Player_Stats shooter_stats)
{
	int i;
	int target_coords[2];
	float max_chance, score_chance, pass_team;

	score_chance = shoot_chance(shooter->final_coords, shooter_stats.post, shooter_stats.shooting);
	max_chance = score_chance;
	target_coords[0] = shooter_stats.post[0] ;
	target_coords[1] = shooter_stats.post[1] ;


	for (i=0; i<shooter_stats.court_group_size-1; i++)
	{
		if (other_players[i].team==shooter->team)
		{
			if (shooter->final_coords[0]!=other_players[i].final_coords[0] && shooter->final_coords[1]!=other_players[i].final_coords[1]) {

				// if my teammate is at a better position than me
				if (score_chance < shoot_chance(other_players[i].final_coords, shooter_stats.post, shooter_stats.shooting)) {
					pass_team = shoot_chance(shooter->final_coords, other_players[i].final_coords, shooter_stats.shooting);
					if (max_chance < pass_team && score_chance < PASS_FACTOR * pass_team) {

						target_coords[0] = other_players[i].final_coords[0] ;
						target_coords[1] = other_players[i].final_coords[1] ;
						max_chance = pass_team;
					}
				}
			}
		}
	}

	// if ball lands 
	if (max_chance >= rand()%100) {
		shooter->ball_coords[0] = target_coords[0] ;
		shooter->ball_coords[1] = target_coords[1] ;
	}

	// negative to indicate that ball did not land
	else {
		shooter->ball_coords[0] = -target_coords[0] ;
		shooter->ball_coords[1] = -target_coords[1] ;
	}
}


void choose_ballpos(int coords[])
{
	int x, y, new_coords[2];  


	if (coords[0] < 0 || coords[1] < 0) {
		// make sure that the randomized ball coods is not the same as target pos
		do {

			x = rand() % 9 ;
			y = rand() % (9-x) ;	

			if (rand()%2==0)	// randomize the direction 
				x = -x;
			if (rand()%2==0)
				y = -y;

			new_coords[0] = -coords[0] + x;
			new_coords[1] = -coords[1] + y;

			limit_coords(new_coords) ;

		} while(abs(coords[0]) == new_coords[0] && abs(coords[1]) == new_coords[1]); 


		coords[0] = new_coords[0] ;
		coords[1] = new_coords[1] ;
	}
}

int main(int argc, char *argv[])
{
	int i,j, winner, field1_winner, numtasks, rounds=1, field_number, game_set;
	int score[3] = {0} ;
	time_t start,end;
	long long before, after ;
	double dif;


	Player_Match match;
	Player_Stats stats;
	MPI_Comm field_comm;
	Player_Match match_info[NPROCS-2];	// -2 field process
	MPI_Datatype MPI_MATCH ;	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

	if (numtasks!=NPROCS) {
		printf("Must specify %d processors. Terminating.\n",NPROCS);
		MPI_Finalize();
		return 1;	
	}


	MPI_MATCH = createPlayerMatchStruct() ;
	MPI_Type_commit(&MPI_MATCH); 

	init_player_stats(&stats);

	init_team(&stats, &match);

	//Timing Code
	/* 
	time(&start);
	before = wall_clock_time() ;
	*/

	for (j=1; j<=2; j++) {
		init_player_match(&match) ;
		for (rounds=1, winner = -1, field1_winner=-1; rounds<=HALF_MATCH; rounds++, winner = -1, field1_winner=-1) {
			// inform players of ball coords
			if (stats.global_rank==FIELD) {
				for (i=1; i<NPROCS; i++)
					MPI_Send(match.ball_coords, 2, MPI_INT, i, T_inform_ball_coords, MPI_COMM_WORLD);
			}
			else {
				MPI_Recv(match.ball_coords, 2, MPI_INT, FIELD, T_inform_ball_coords, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

			if (stats.global_rank!=FIELD_0 && stats.global_rank!=FIELD_1) {
				player_run(&match, &stats) ;
			}



			setup_field_comm(&stats, &match, &field_comm); 

			// send coordinates together with other info
			if (stats.court_rank==FIELD) {
				for (i=1; i<stats.court_group_size; i++)
					MPI_Recv(&match_info[i-1], 1, MPI_MATCH, i, T_report_pos, field_comm, MPI_STATUS_IGNORE);
			}
			else {
				MPI_Send(&match, 1, MPI_MATCH, 0, T_report_pos, field_comm);
			}

			// decide winner and inform players. Players in other field will not know who is the winner
			if (stats.court_rank==FIELD) {

				winner = decide_winner(match_info, stats.court_group_size-1) ;
				for (i=1; i<stats.court_group_size; i++)
				{
					MPI_Send(&winner, 1, MPI_INT, i, T_inform_winner, field_comm);

				}
			}
			else
			{
				MPI_Recv(&winner, 1, MPI_INT, FIELD, T_inform_winner, field_comm, MPI_STATUS_IGNORE);

			}


			if (winner!=-1) {

				if (stats.court_rank==winner) {
					// winner receives information of all players in same field 
					MPI_Recv(match_info, stats.court_group_size-1, MPI_MATCH, FIELD, T_players_on_field, field_comm, MPI_STATUS_IGNORE);
					shoot(match_info, &match, stats);

					// winner sends shoot location (negative for unsucessful shoot)
					MPI_Send(match.ball_coords, 2, MPI_INT, 0, T_inform_shoot_pos, field_comm);
				}

				if (stats.court_rank==FIELD)
				{	
					MPI_Send(match_info, stats.court_group_size-1, MPI_MATCH, winner, T_players_on_field, field_comm);
					MPI_Recv(match.ball_coords, 2, MPI_INT, winner, T_inform_shoot_pos, field_comm, MPI_STATUS_IGNORE);

					match_info[winner-1].ball_coords[0] = match.ball_coords[0] ;
					match_info[winner-1].ball_coords[1] = match.ball_coords[1] ;

					// decides the final position of ball
					choose_ballpos(match.ball_coords) ;
				}


			}

			// field 1 send to field 0 if it has more than 1 player in its court (not including itself)
			if (stats.global_rank==FIELD_1 && stats.court_group_size > 1) {
				MPI_Send(match_info, stats.court_group_size - 1, MPI_MATCH, FIELD_0, T_update_field0_player, stats.team_comm);

			}
			//field 0 receive from field 1 if it does not have all the players in its court (total processes - process of field 1)
			if (stats.global_rank==FIELD_0 && stats.court_group_size != NPROCS-1) {
				// buffer has an offset of no. of players in its court so that its own players will not be overwritten
				MPI_Recv(&match_info[stats.court_group_size - 1], NPROCS - stats.court_group_size -1, MPI_MATCH, FIELD_1, T_update_field0_player, stats.team_comm, MPI_STATUS_IGNORE);
			}

			// send the winner of field 1 to field 0
			if (stats.global_rank==FIELD_1) {
				MPI_Send(&winner, 1, MPI_INT, FIELD_0, T_winner_field1, MPI_COMM_WORLD);
			}
			if (stats.global_rank==FIELD_0) { 
				MPI_Recv(&field1_winner, 1, MPI_INT, FIELD_1, T_winner_field1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

			// field 0 update ball coords there is a winner from field 1 
			if (winner!=-1 && stats.global_rank==FIELD_1) {
				MPI_Send(match.ball_coords, 2, MPI_INT, FIELD_0, T_update_field0_ball, MPI_COMM_WORLD);
			}
			if (winner==-1 && stats.global_rank==FIELD_0 && field1_winner!=-1)
				MPI_Recv(match.ball_coords, 2, MPI_INT, FIELD_1, T_update_field0_ball, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			if (stats.global_rank==FIELD_0) {
				printf("%d\n", rounds*j) ; 
				print_round_info(match_info, score, match.ball_coords);

			}
			else {
				match.ball_coords[0] = INVALID_COORDS;
				match.ball_coords[1] = INVALID_COORDS;
				match.initial_coords[0] = match.final_coords[0] ;
				match.initial_coords[1] = match.final_coords[1] ;
			}

			MPI_Comm_free(&field_comm);
		}

		if (match.team == 1) {
			stats.post[0] = 128;
		}
		else if  (match.team == 2) {
			stats.post[0] = 0;
		}
	}

	//Timing Code
	/*
	if (stats.global_rank==FIELD_0) {
	time(&end);	
	/dif = difftime (end,start);
	after = wall_clock_time() ;
	printf ("It took %.5lf seconds to complete 5400 rounds (time.h).\n", dif );
	printf ("It took %.5lf seconds to complete 5400 rounds (linux time).\n", (float)(after-before)/1000000000.0);
	printf ("It took %.5lf milliseconds to complete 1 rounds (linux time).\n", (float)(after-before)/(5400*1000000.0));
	}
	*/
	MPI_Finalize();
	return 0;
}
