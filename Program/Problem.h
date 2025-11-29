#ifndef _PROBLEM_H
#define _PROBLEM_H

#define INFINITO 999999999

#include <math.h>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "Data.h"

//------ DEFINITION OF TYPES OF PROBLEM SPECIFIC --------

struct TDelivery
{
 int value;
 int origin;
 int destination;
};

//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  --------

static std::vector <std::vector <int> > dist;	// matrix with Euclidean distance
static std::vector <TDelivery> deliveries;
static int capacity;



//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------


/************************************************************************************
 Method: ReadData
 Description: read input data of the problem
*************************************************************************************/
void ReadData(char nameTable[], int &n);

/************************************************************************************
 Method: Decoder()
 Description: Convert a random key solution in a real problem solution
*************************************************************************************/
void Decoder(TSol &s, int n, int nDec);

/************************************************************************************
 Method: LocalSearch
 Description: RVND
*************************************************************************************/
void LocalSearch(TSol &s, int n, int nLS);

/************************************************************************************
 Method: LocalSearchExactDeliveries
 Description: RVND
*************************************************************************************/
void LocalSearchExactDeliveries(TSol &s, int n, int nLS);


/************************************************************************************
 Method: CalculateFitness
 Description: calculate the fitness of a chromossom s
*************************************************************************************/
double CalculateFitness(TSol &s, int n, long nodeSize);

/************************************************************************************
 Method: Dec1
 Description: sort decoder 
*************************************************************************************/
void Dec1(TSol &s, int n);

/************************************************************************************
 Method: Dec2
 Description: 2-opt decoder 
*************************************************************************************/
void Dec2(TSol &s, int n);

/************************************************************************************
 Method: Dec3
 Description: Cheapest Insertion decoder 
*************************************************************************************/
void Dec3(TSol &s, int n);

/************************************************************************************
 Method: Dec4
 Description: k-Farthest Insertion decoder 
*************************************************************************************/
void Dec4(TSol &s, int n);

/************************************************************************************
 Method: Dec5
 Description: k-Nearest Insertion decoder 
*************************************************************************************/
void Dec5(TSol &s, int n);

/************************************************************************************
 Method: LS1
 Description: NodeExchange
*************************************************************************************/
void LS1(TSol &s, std::vector<int> &sortedDeliveriesByValue);

/************************************************************************************
 Method: LS1Exact
 Description: NodeExchange
*************************************************************************************/
void LS1Exact(TSol &s);

/************************************************************************************
 Method: LS2
 Description: NodeInsertion
*************************************************************************************/
void LS2(TSol &s, std::vector<int> &sortedDeliveriesByValue);

/************************************************************************************
 Method: LS2Exact
 Description: NodeInsertion
*************************************************************************************/
void LS2Exact(TSol &s);

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem();

#endif

void printSolution(TSol &s, int n);

void sanityCheck(TSol s);