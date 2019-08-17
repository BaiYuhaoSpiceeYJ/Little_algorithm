/**
 * Projec : gtsp (voyageur de commerce)
 *
 * Date   : 07/04/2014
 * Author : Olivier Grunder
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NBR_TOWNS 6
/* Distance matrix */
double dist[NBR_TOWNS][NBR_TOWNS] ;

/* Each edge has a starting and ending node */
int starting_town[NBR_TOWNS] ;
int ending_town[NBR_TOWNS] ;

/* no comment */
int best_solution[NBR_TOWNS] ;
double best_eval=-1.0 ;
int nbr_operation = 0;

/**
 * Berlin52 :
 *  6 towns : Best solution (2315.15): 0 1 2 3 5 4
 * 10 towns : Best solution (2826.50): 0 1 6 2 7 8 9 3 5 4
 */
float coord[NBR_TOWNS][2]= {
	{565.0,  575.0},
	{ 25.0,  185.0},
	{345.0,  750.0},
	{945.0,  685.0},
	{845.0,  655.0},
	{880.0,  660.0},
	//{ 25.0,  230.0},
	//{525.0,  1000.0},
	//{580.0,  1175.0},
	//{650.0,  1130.0},
} ;



/**
 * print a matrix
 */
void print_matrix(double d[NBR_TOWNS][NBR_TOWNS]) {
	int i, j ;
	for (i=0; i<NBR_TOWNS; i++) {
		printf ("%d:", i) ;
		for (j=0; j<NBR_TOWNS; j++) {
			printf ("%6.1f ", d[i][j]) ;
		}
		printf ("\n") ;
	}
}



/**
 * print a solution
 */
void print_solution(int* sol, double eval) {
	int i ;
	printf ("(%.2f): ", eval) ;
	for (i=0; i<NBR_TOWNS; i++)
		printf ("%d ",sol[i]);
	printf("\n") ;
}




/**
 * evaluation of a solution
 */
double evaluation_solution(int* sol) {
	double eval=0 ;
	int i ;
	for (i=0; i<NBR_TOWNS-1; i++) {
		eval += dist[sol[i]][sol[i+1]] ;
	}
	eval += dist[sol[NBR_TOWNS-1]][sol[0]] ;

	return eval ;

}




/**
 * nearest neighbour solution
 */
double build_nearest_neighbour() {
	/* solution of the nearest neighbour */
	int i,j, sol[NBR_TOWNS] ;

	/* evaluation of the solution */
	double eval = 0 ;

	sol[0] = 0 ;

	int visited[NBR_TOWNS];
	for (i=0; i<NBR_TOWNS; i++) visited[i] = 0;
	visited[0] = 1;
	int cur_town = 0;


	for (i=0; i<NBR_TOWNS-1; i++) {
		double nearest_dist = 9999.9;
		int nearest_neighbor = cur_town;
		for (j=0; j<NBR_TOWNS; j++) {
			if(visited[j] == 0 && dist[cur_town][j] < nearest_dist && j!=cur_town) {
				nearest_neighbor = j;
				nearest_dist = dist[cur_town][j];
			}
		}
		visited[nearest_neighbor] = 1;
		sol[i+1] = nearest_neighbor;
		cur_town = nearest_neighbor;
	}

	eval = evaluation_solution(sol) ;
	printf("Nearest neighbour ") ;
	print_solution (sol, eval) ;

	for (i=0; i<NBR_TOWNS; i++) best_solution[i] = sol[i] ;
	best_eval  = eval ;


	return eval ;
}




/**
 *  Build final solution
 */
void build_solution() {
	int i, solution[NBR_TOWNS] ;

	int indiceCour = 0;
	int villeCour = 0;

	while (indiceCour < NBR_TOWNS) {

		solution[indiceCour] = villeCour ;

		// Test si le cycle est hamiltonien
		for (i = 0; i < indiceCour; i++) {
			if (solution[i] == villeCour) {
				/* printf ("cycle non hamiltonien\n") ; */
				return;
			}
		}
		// Recherche de la ville suivante
		int trouve = 0;
		int i = 0;
		while ((!trouve) && (i < NBR_TOWNS)) {
			if (starting_town[i] == villeCour) {
				trouve = 1;
				villeCour = ending_town[i];
			}
			i++;
		}
		indiceCour++;
	}

	double eval = evaluation_solution(solution) ;

	if (best_eval<0 || eval < best_eval) {
		best_eval = eval ;
		for (i=0; i<NBR_TOWNS; i++)
			best_solution[i] = solution[i] ;
		printf ("New best solution: ") ;
		print_solution (solution, best_eval) ;
	}
	return;
}




/**
 *  Little Algorithm
 */
void little_algorithm(double d0[NBR_TOWNS][NBR_TOWNS], int iteration, double eval_node_parent) {
	nbr_operation++;
	printf("best eval = %f\n",best_eval);
	printf("iteration%d\n",iteration);
	if (iteration >= NBR_TOWNS+1) {
		build_solution ();
		return;
	}

	/* Do the modification on a copy of the distance matrix */
	double d[NBR_TOWNS][NBR_TOWNS] ;
	memcpy (d, d0, NBR_TOWNS*NBR_TOWNS*sizeof(double)) ;

	int i, j ;

	double eval_node_left_child = eval_node_parent;
	double eval_node_right_child = eval_node_parent;

	/**
	 * substract the min of the rows and the min of the columns
	 * and update the evaluation of the current node
	 */
	printf("Substract the min of the rows begin\n");
	int total_substract_line = 0;
	for(i=0; i<NBR_TOWNS; i++) {
		double min = 9999.9;
		for(j=0; j<NBR_TOWNS; j++) {
			if(d[i][j] < min && i != j && d[i][j]!=9999.9) {
				min = d[i][j];
			}
		}
		total_substract_line += min;
		for(j=0; j<NBR_TOWNS; j++) {
			if(d[i][j]!=9999.9) d[i][j] -= min;
		}

	}
	print_matrix (d);
	printf("Substract the min of the rows finished\n\n");

	printf("Substract the min of the cols begin\n");
	int total_substract_col = 0;
	for(i=0; i<NBR_TOWNS; i++) {
		double min = 9999.9;
		for(j=0; j<NBR_TOWNS; j++) {
			if(d[j][i] < min && i != j && d[j][i]!=9999.9) {
				min = d[j][i];
			}
		}
		total_substract_col += min;
		for(j=0; j<NBR_TOWNS; j++) {
			if(d[j][i]!=9999.9)d[j][i] -= min;
		}
	}
	eval_node_left_child += total_substract_col + total_substract_line;
	print_matrix (d);
	printf("Substract the min of the cols finished\n\n");


	/* Cut : stop the exploration of this node */
	if (best_eval>=0 && eval_node_left_child >= best_eval)
		return;


	/**
	 *  Compute the penalities to identify the zero with max penalty
	 *  If no zero in the matrix, then return, solution infeasible
	/* row and column of the zero with the max penalty */
	double max_penalty = 0;
	int max_penalty_line = -1;
	int max_penalty_col = -1;
	for(i=0; i<NBR_TOWNS; i++) {
		for(j=0; j<NBR_TOWNS; j++) {
			double cur_penalty_line = 9999.9;
			double cur_penalty_col = 9999.9;
			double cur_penalty;
			int k;
			if(d[i][j]!=0) continue;
			else {
				for(k=0; k<NBR_TOWNS; k++) {
					if(d[i][k]<cur_penalty_line && d[i][k]!=9999.9 && j!=k) {
						cur_penalty_line = d[i][k];
					}
					if(d[k][j]<cur_penalty_col && d[k][j]!=9999.9 && i!=k) {
						cur_penalty_col = d[k][j];
					}
				}
			}
			cur_penalty = cur_penalty_line + cur_penalty_col;
			printf("current penalty=%f ,current row=%d ,current col=%d\n",cur_penalty,i,j);
			if(cur_penalty >= max_penalty) {
				max_penalty = cur_penalty;
				max_penalty_line = i;
				max_penalty_col = j;
			}
		}
	}
	if(max_penalty_line==-1||max_penalty_col==-1) return;
	printf("max penalty=%f ,max_line=%d ,max_col=%d\n\n",max_penalty,max_penalty_line,max_penalty_col);


	/**
	 *  Store the row and column of the zero with max penalty in
	 *  starting_town and ending_town
	 */
	int izero = max_penalty_line;//starting_town
	int jzero = max_penalty_col;//ending_town
	starting_town[iteration-1] = izero;
	ending_town[iteration-1] = jzero;
	eval_node_right_child += max_penalty;
	/* Do the modification on a copy of the distance matrix */
	double d2[NBR_TOWNS][NBR_TOWNS] ;
	memcpy (d2, d, NBR_TOWNS*NBR_TOWNS*sizeof(double)) ;

	/**
	 *  Modify the matrix d2 according to the choice of the zero with the max penalty
	 */
	d2[izero][jzero] = 9999.9;
	d2[jzero][izero] = 9999.9;
	for(i=0; i<NBR_TOWNS; i++) {
		for(j=0; j<NBR_TOWNS; j++) {
			if(i==max_penalty_line || j==max_penalty_col) {
				d2[i][j] = 9999.9;
			} else {
				d2[i][j] = d[i][j];
			}

		}
	}
	printf("The new matrice is\n");
	print_matrix (d2);
	printf("\niteration %d end\n",iteration);
	for(i=0; i<NBR_TOWNS; i++) {
		printf("for this iteration, the starting town is %d and the ending town is %d\n",starting_town[i],ending_town[i]);
	}
	printf("\n\n\n\n\n\n");

	/* Explore left child node according to given choice */
	little_algorithm(d2, iteration+1 , eval_node_left_child);

	/* Do the modification on a copy of the distance matrix */
	memcpy (d2, d, NBR_TOWNS*NBR_TOWNS*sizeof(double)) ;
	printf("\n\n\n\n\n\nchoose the non-choice of %d-%d\n",izero,jzero);


	/**
	 *  Modify the dist matrix to explore the other possibility : the non-choice
	 *  of the zero with the max penalty
	Explore right child node according to non-choice */
	d2[izero][jzero] = 9999.9;
	d2[jzero][izero] = 9999.9;
	printf("Substract the min of the rows begin\n");
	for(i=0; i<NBR_TOWNS; i++) {
		double min = 9999.9;
		for(j=0; j<NBR_TOWNS; j++) {
			if(d2[i][j] < min && i != j && d2[i][j]!=9999.9) {
				min = d2[i][j];
			}
		}
		for(j=0; j<NBR_TOWNS; j++) {
			if(d2[i][j]!=9999.9) d2[i][j] -= min;
		}

	}
	print_matrix (d2);
	printf("Substract the min of the rows finished\n\n");

	printf("Substract the min of the cols begin\n");
	for(i=0; i<NBR_TOWNS; i++) {
		double min = 9999.9;
		for(j=0; j<NBR_TOWNS; j++) {
			if(d2[j][i] < min && i != j && d2[j][i]!=9999.9) {
				min = d2[j][i];
			}
		}
		for(j=0; j<NBR_TOWNS; j++) {
			if(d2[j][i]!=9999.9)d2[j][i] -= min;
		}
	}
	print_matrix (d2);
	//return;
	printf("Substract the min of the cols finished\n\n");
	if (best_eval>=0 && eval_node_left_child >= best_eval)
		return;
	little_algorithm(d2, iteration, eval_node_right_child);

}




/**
 *
 */
int main (int argc, char* argv[]) {

	best_eval = -1 ;
	/* Print problem informations */
	printf ("Points coordinates:\n") ;
	int i ;
	int j;
	for (i=0; i<NBR_TOWNS; i++) {
		printf ("Node %d: x=%f, y=%f\n", i, coord[i][0], coord[i][1]) ;
	}
	printf ("\n") ;


	/* Calcul de la matrice des distances */
	for(i=0; i<NBR_TOWNS; i++) {
		for(j=0; j<NBR_TOWNS; j++) {
			dist [i][j] = sqrt((coord[i][0]-coord[j][0]) * (coord[i][0]-coord[j][0])+ (coord[i][1]-coord[j][1]) * (coord[i][1]-coord[j][1]));
		}
		dist[i][i] = 9999.9;
	}

	printf ("Distance Matrix:\n") ;
	print_matrix (dist) ;
	printf ("\n") ;


	double nearest_neighbour = build_nearest_neighbour() ;
	printf("%6.1f ", nearest_neighbour);
	printf("best eval = %f\n",best_eval);
	best_eval = -1 ;
	little_algorithm(dist,1,0);


	printf("Best solution:") ;
	print_solution (best_solution, best_eval) ;
	printf("total number of operation = %d\n",nbr_operation);
	printf ("Hit RETURN!\n") ;
	getchar() ;


	return 0 ;
}
