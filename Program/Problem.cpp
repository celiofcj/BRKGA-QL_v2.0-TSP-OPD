#include "Problem.h"

#include <map>

#include "LazyPropagationSegmentTree.h"
#include "MCMF.h"
#include "Searcher.h"


// Sort TSol by random-keys
bool sortByRk(const TVecSol &lhs, const TVecSol &rhs) { return lhs.rk < rhs.rk; }

void ReadData(char nameTable[], int &n)
{
    capacity = 10;
    char name[200] = "../Instances/";
    strcat(name,nameTable);

    FILE *arq;
    arq = fopen(name,"r");

    if (arq == NULL)
    {
        printf("\nERROR: File (%s) not found!\n",name);
        getchar();
        exit(1);
    }

    char line[256];
    int nodesSize = 0, deliveriesSize = 0;

    fgets(line, sizeof(line), arq);
    sscanf(line,"%d %d",&nodesSize, &deliveriesSize);

    n = nodesSize + deliveriesSize;

    dist.clear();
    dist.resize(nodesSize, std::vector<int>(nodesSize));

    for (int i = 0; i < nodesSize; i++)
    {
        fgets(line, sizeof(line), arq);
        char *p = line;
        int temp = -1;
        int j = 0;
        while (sscanf(p, "%d", &temp) == 1) {
            if (j == nodesSize) break;
            dist[i][j] = temp;
            j++;

            while (*p && *p != ' ') p++;
            p++;
        }
    }

    TDelivery deliveryTemp;
    for (int i = 0; i < deliveriesSize; i++)
    {
        fgets(line, sizeof(line), arq);
        sscanf(line, "%d %d %d", &deliveryTemp.origin, &deliveryTemp.destination, &deliveryTemp.value);
        deliveries.push_back(deliveryTemp);
    }

    fclose(arq);

    // for (std::vector<int> node : dist)
    // {
    //     for (int i : node)
    //     {
    //         printf("%d ", i);
    //     }
    //     printf("\n");
    // }
    //
    // for (TDelivery delivery : deliveries) {
    //     printf("%d %d %d\n", delivery.origin, delivery.destination, delivery.value);
    // }
}

void Decoder(TSol &s, int n, int nDec)
{
    // copy the random-key sequence of current solution
    TSol temp = s;

    // define decoder function based in the random-key of position n+1
    int dec = floor(s.vec[n].rk*nDec)+1;
    switch (dec)
    {
        case 1:
            Dec1(s, n);
            break;

        case 2:
            Dec2(s, n);
            break;

        default:
            break;
    }

    // return initial random-key sequence and maintain the solution sequence
    for (int i=0; i<n; i++){
        s.vec[i].rk = temp.vec[i].rk;
    }
}


bool sortByIncreasingValue(const std::pair<int, TDelivery>& a, const std::pair<int, TDelivery>& b) {
    return a.second.value > b.second.value;
}

std::vector<int> ObtainSortedDeliveries(TSol &s, int n, int routeLength) {
    int deliverySize = n - routeLength;

    std::map<int, TDelivery> map;
    std::vector<int> newVector(deliverySize);

    for (int i = 0; i < deliverySize; i++) {
        int index = s.vec[routeLength + i].sol;

        if (index < 0) {
            index = -(index + 1);
        }


        map[index] = deliveries[index];
    }

    std::vector<std::pair<int, TDelivery>> vec(map.begin(), map.end());

    sort(vec.begin(), vec.end(), sortByIncreasingValue);

    int i = 0;
    for (const auto& pair : vec) {
        newVector[i] = pair.first;
        i++;
    }

    return newVector;
}


void LocalSearch(TSol &s, int n, int nLS)
{
    // ***** we use a Random Variable Neighborhood Descent (RVND) as local search ****

    int routeLength = dist.size();
    std::vector<int> sortedDeliveriesByValue = ObtainSortedDeliveries(s, n, routeLength);

	int k = 1;
    std::vector <int> NSL;
    std::vector <int> NSLAux;

    for (int i=1; i<=nLS; i++)
    {
        NSL.push_back(i);
        NSLAux.push_back(i);
    }

	while (!NSL.empty())
	{
        // current objective function
        double foCurrent = s.ofv;

        // randomly choose a neighborhood
        int pos = rand() % NSL.size();
        k = NSL[pos];

        switch (k)
        {
            case 1:
                LS1(s, sortedDeliveriesByValue);
                break;

            case 2:
                LS2(s, sortedDeliveriesByValue);
                break;
            default:
                break;
        }

        // we better the current solution
        if (s.ofv < foCurrent)
        {
            // refresh NSL
            NSL.clear();
            NSL = NSLAux;
        }
        else
        {
            // Remove N(n) from NSL
            NSL.erase(NSL.begin()+pos);
        }
    } //end while

    for (int i=1; i<=nLS; i++)
    {
        NSL.push_back(i);
        NSLAux.push_back(i);
    }

    while (!NSL.empty())
    {
        // current objective function
        double foCurrent = s.ofv;

        // randomly choose a neighborhood
        int pos = rand() % NSL.size();
        k = NSL[pos];

        switch (k)
        {
            case 1:
                LS1Exact(s);
                break;

            case 2:
                LS2Exact(s);
                break;
            default:
                break;
        }

        // we better the current solution
        if (s.ofv < foCurrent)
        {
            // refresh NSL
            NSL.clear();
            NSL = NSLAux;
        }
        else
        {
            // Remove N(n) from NSL
            NSL.erase(NSL.begin()+pos);
        }
    } //end while


}

double CalculateFitness(TSol &s, int n, long nodeSize)
{
    s.ofv = 0;
    s.ofvD = 0;
    s.ofvR = 0;

    for (int i=0 ; i < nodeSize; i++)
     {
        s.ofvR += dist[s.vec[i%nodeSize].sol][s.vec[(i+1)%nodeSize].sol];
     }

    for (int i=nodeSize; i < n; i++)
    {
        if (s.vec[i].sol >= 0) {
            s.ofvD += deliveries[s.vec[i].sol].value;
        }
    }

    s.ofv = s.ofvR - s.ofvD;

    return s.ofv;
}

double SumDeliveries(std::vector<int> del) {
    double sum = 0;
    for (int i = 0; i < del.size(); i++) {
        if (del[i] >= 0) {
            sum += deliveries[i].value;
        }
    }

    return sum;
}

void InitSolutionNodes(TSol &s, long size) {
    s.vec[0].sol = 0;
    for (int i = 1; i < size; i++)
    {
        s.vec[i].sol = i;
    }
}

void InitSolutionDeliveries(TSol &s, int n, long startingIndex) {
    for (int i = 0; i + startingIndex < n; i++) {
        s.vec[startingIndex + i].sol = i;
    }
}

bool withinCapacity(int origin, int destination, int nodesQtd, LazyPropagationSegmentTree &segment_tree) {
    if (destination == 0)
    {
        destination = nodesQtd;
    }

    return segment_tree.query(origin, destination - 1) < capacity;
}

void updateWeight(int origin, int destination, int nodesQtd, int addend, LazyPropagationSegmentTree &segment_tree) {
    if (destination == 0)
    {
        destination = nodesQtd;
    }

    segment_tree.update(origin, destination - 1, addend);
}

Searcher GetSearcher(TSol &s, long nodesSize) {
    std::vector<int> nodeSols = std::vector<int>(nodesSize);
    for (int i = 0; i < nodesSize; i++)
    {
        nodeSols[i] = s.vec[i].sol;
    }

    return {nodesSize, nodeSols};
}

bool VerifyDeliveryConstrains(TSol &s, long routeLength, LazyPropagationSegmentTree &segment_tree, Searcher searcher, int it) {
    TDelivery delivery = deliveries[s.vec[it].sol];

    // printf("Delivery %ld: id: %d - pos: %d - value: %d - origin: %d - destination: %d\n", i-nodesSize, s.vec[i].sol, s.vec[i].sol, delivery.value, delivery.origin, delivery.destination);

    int origin = searcher.search(delivery.origin);
    int destination = searcher.search(delivery.destination);

    // printf("I: %d, Origin: %d, Destination: %d\n", i, origin, destination);

    if ((origin < destination || (destination == 0 && origin <= routeLength - 1)) && withinCapacity(origin, destination, routeLength, segment_tree))
    {
        updateWeight(origin, destination, routeLength, 1, segment_tree);
        return true;
    }

    return false;
}

bool VerifyDeliveryConstrainsRaw(int routeLength, std::vector<int> &del, LazyPropagationSegmentTree &segment_tree, Searcher &searcher, int &it) {
    TDelivery delivery = deliveries[del[it]];

    // printf("Delivery %ld: id: %d - pos: %d - value: %d - origin: %d - destination: %d\n", i-nodesSize, s.vec[i].sol, s.vec[i].sol, delivery.value, delivery.origin, delivery.destination);

    int origin = searcher.search(delivery.origin);
    int destination = searcher.search(delivery.destination);

    // printf("I: %d, Origin: %d, Destination: %d\n", i, origin, destination);

    if ((origin < destination || (destination == 0 && origin <= routeLength - 1)) && withinCapacity(origin, destination, routeLength, segment_tree))
    {
        updateWeight(origin, destination, routeLength, 1, segment_tree);
        return true;

    }

    return false;
}

void SelectDeliveries(TSol &s, int n, long nodesSize) {
    LazyPropagationSegmentTree segment_tree = LazyPropagationSegmentTree(nodesSize + 1);
    Searcher searcher = GetSearcher(s, nodesSize);

    for (int i = nodesSize; i < n; i++)
    {
        if (!VerifyDeliveryConstrains(s, nodesSize, segment_tree, searcher, i)) {
            s.vec[i].sol = -s.vec[i].sol - 1;
        }
    }
}

int EvaluateDeliveries(TSol &s, long routeLength, std::vector<int> &del) {
    LazyPropagationSegmentTree segment_tree = LazyPropagationSegmentTree(routeLength + 1);
    Searcher searcher = GetSearcher(s, routeLength);

    int sum = 0;
    for (int i = 0; i < del.size(); i++)
    {
        if (VerifyDeliveryConstrainsRaw(routeLength, del, segment_tree, searcher, i)) {
            sum += deliveries[del[i]].value;
        }
    }

    return sum;
}

void Dec1(TSol &s, int n) // sort
{
    long nodesSize = dist.size();

    InitSolutionNodes(s, nodesSize);

    InitSolutionDeliveries(s, n, nodesSize);

    sort(s.vec.begin() + 1, s.vec.begin() + nodesSize, sortByRk);

    sort(s.vec.begin() + nodesSize + 1, s.vec.end()-1, sortByRk);

    SelectDeliveries(s, n, nodesSize);

    s.vec[n].sol = -1;

    CalculateFitness(s,n, nodesSize);
}

void SelectDeliveriesBinary(TSol &s, int n, long nodeSize)
{
    LazyPropagationSegmentTree segment_tree = LazyPropagationSegmentTree(nodeSize + 1);
    Searcher searcher = GetSearcher(s, nodeSize);

    int deliverySize = n - nodeSize;
    int startPosition = floor(s.vec[0].rk * (deliverySize));
    for (int i = 0; i < deliverySize; i++)
    {
        int it = (startPosition + i) % deliverySize + nodeSize;
        if (s.vec[it].rk <= 0.5 || !VerifyDeliveryConstrains(s, nodeSize, segment_tree, searcher, it)) {
            s.vec[it].sol = -s.vec[it].sol - 1;
        }
    }
}

void Dec2(TSol &s, int n)
{
    long nodesSize = dist.size();

    InitSolutionNodes(s, nodesSize);

    InitSolutionDeliveries(s, n, nodesSize);

    sort(s.vec.begin() + 1, s.vec.begin() + nodesSize, sortByRk);

    SelectDeliveriesBinary(s, n, nodesSize);

    CalculateFitness(s,n, nodesSize);
}

std::vector<int> OptimumDeliveries(TSol &s) {
    MCMF mcmf(dist.size(), deliveries.size());

    mcmf.build(s, deliveries, capacity);

    int del[deliveries.size()];

    for (int i = 0; i < deliveries.size(); i++) {
        del[i] = -1;
    }

    mcmf.solve(del);

    std::vector<int> optimumDeliveries(deliveries.size());

    for (int i = 0; i < deliveries.size(); i++) {
        if (del[i] == 1) {
            optimumDeliveries[i] = i;
        } else {
            optimumDeliveries[i] = -i - 1;
        }
    }

    return optimumDeliveries;
}

void SwitchNodes(TSol &s, const int &i, const int &j) {
    TVecSol aux = s.vec[i];

    s.vec[i] = s.vec[j];
    s.vec[j] = aux;
}

void LS1(TSol &s, std::vector<int> &sortedDeliveriesByValue) // NodeExchange
{
    int bestI = -1, bestJ = -1;
    double bestFoOptD = 0, bestFoOptR = 0, bestFoOpt = 0;

    int routeLength = dist.size();

    for (int i = 1; i < routeLength - 1; i++)
    {
        for (int j = i + 2; j < routeLength; j++)
        {
            if (j != routeLength - 1) //no exchange edge of the tour
            {
                int vi = s.vec[i].sol;
                int viP = s.vec[i+1].sol;
                int viM = s.vec[i-1].sol;

                int vj = s.vec[j].sol;
                int vjM = s.vec[j-1].sol;
                int vjP = 0;
                if (j < routeLength - 1)
                    vjP = s.vec[j+1].sol;
                else
                    vjP = s.vec[0].sol;

                double foOptR = - dist[viM][vi]
                                - dist[vi][viP]
                                - dist[vjM][vj]
                                - dist[vj][vjP]
                                + dist[viM][vj]
                                + dist[vj][viP]
                                + dist[vjM][vi]
                                + dist[vi][vjP];

                SwitchNodes(s, i, j);

                double foOptD = s.ofvD - EvaluateDeliveries(s, routeLength, sortedDeliveriesByValue);

                SwitchNodes(s, i, j);

                double foOpt = foOptR + foOptD;

                if (foOpt < bestFoOpt) {
                    bestFoOpt = foOpt;
                    bestFoOptR = foOptR;
                    bestFoOptD = foOptD;
                    bestI = i;
                    bestJ = j;
                }
            }
        }
    }

    double ofv = s.ofv + bestFoOpt;

    if (bestI != -1 && bestJ != -1 && ofv < s.ofv)
    {
        SwitchNodes(s, bestI, bestJ);

        int deliveriesSize = deliveries.size();
        for (int i = 0; i < deliveriesSize; i++)
        {
            s.vec[i + routeLength].sol = sortedDeliveriesByValue[i];
        }

        s.ofv = ofv;
        s.ofvR += bestFoOptR;
        s.ofvD -= bestFoOptD;

        // sanityCheck(s);
    }
}

void LS1Exact(TSol &s) {
    int bestI = -1, bestJ = -1;
    double bestFoOptD = 0, bestFoOptR = 0, bestFoOpt = 0;
    std::vector<int> bestDeliveries(deliveries.size());

    int routeLength = dist.size();

    for (int i = 1; i < routeLength - 1; i++)
    {
        for (int j = i + 2; j < routeLength; j++)
        {
            if (j != routeLength - 1) //no exchange edge of the tour
            {
                int vi = s.vec[i].sol;
                int viP = s.vec[i+1].sol;
                int viM = s.vec[i-1].sol;

                int vj = s.vec[j].sol;
                int vjM = s.vec[j-1].sol;
                int vjP = 0;
                if (j < routeLength - 1)
                    vjP = s.vec[j+1].sol;
                else
                    vjP = s.vec[0].sol;

                double foOptR = - dist[viM][vi]
                                - dist[vi][viP]
                                - dist[vjM][vj]
                                - dist[vj][vjP]
                                + dist[viM][vj]
                                + dist[vj][viP]
                                + dist[vjM][vi]
                                + dist[vi][vjP];

                SwitchNodes(s, i, j);

                std::vector<int> optimumDeliveries = OptimumDeliveries(s);

                double foOptD = s.ofvD - SumDeliveries(optimumDeliveries);

                SwitchNodes(s, i, j);

                double foOpt = foOptR + foOptD;

                if (foOpt < bestFoOpt) {
                    bestFoOpt = foOpt;
                    bestFoOptR = foOptR;
                    bestFoOptD = foOptD;
                    bestI = i;
                    bestJ = j;
                    bestDeliveries = optimumDeliveries;
                }
            }
        }
    }

    double ofv = s.ofv + bestFoOpt;

    if (bestI != -1 && bestJ != -1 && ofv < s.ofv)
    {
        SwitchNodes(s, bestI, bestJ);

        int deliveriesSize = deliveries.size();
        for (int i = 0; i < deliveriesSize; i++)
        {
            s.vec[i + routeLength].sol = bestDeliveries[i];
        }

        s.ofv = ofv;
        s.ofvR += bestFoOptR;
        s.ofvD -= bestFoOptD;

        // sanityCheck(s);
    }
}

TSol InsertNode(TSol s, int n, int i, int j) {
    TVecSol aux;

    aux.sol = s.vec[i].sol;
    aux.rk = s.vec[i].rk;

    if (i < j)
    {
        s.vec.insert(s.vec.begin() + j + 1, aux);
        s.vec.erase(s.vec.begin() + i);
    }
    else
    {
        s.vec.insert(s.vec.begin() + j, aux);
        s.vec.erase(s.vec.begin() + i + 1);
    }


    return TSol(s);
}

void LS2(TSol &s,  std::vector<int> &sortedDeliveriesByValue) // NodeInsertion
{
    int bestI = -1, bestJ = -1;
    double bestFoOptD = 0, bestFoOptR = 0, bestFoOpt = 0;

    int routeLength = dist.size();

    for (int i = 1; i < routeLength; i++)
    {
        for (int j = (i + 1) % routeLength; j != (i + 1) % routeLength - 1 && j != i; j = (j + 1) % routeLength)
        {
            if (j == 0) {
                continue;
            }

            int vi  = s.vec[i].sol;
            int viP = s.vec[(i+1)%routeLength].sol;
            int viM = s.vec[i-1].sol;

            int vjM = s.vec[j - 1].sol;
            int vj  = s.vec[j].sol;
            int vjP = s.vec[(j+1)%routeLength].sol;

            double foOptR;

            if (i < j) {
                foOptR = - dist[vj][vjP]
                    - dist[viM][vi]
                    - dist[vi][viP]
                    + dist[vj][vi]
                    + dist[vi][vjP]
                    + dist[viM][viP];
            }
            else
            {
                foOptR = -dist[vjM][vj]
                    - dist[viM][vi]
                    - dist[vi][viP]
                    + dist[vi][vj]
                    + dist[vjM][vi]
                    + dist[viM][viP];

            }

            TSol newSol = InsertNode(s, routeLength, i, j);

            double foOptD = s.ofvD - EvaluateDeliveries(newSol, routeLength, sortedDeliveriesByValue);

            double foOpt = foOptR + foOptD;

            if (foOpt < bestFoOpt) {
                bestFoOpt = foOpt;
                bestFoOptR = foOptR;
                bestFoOptD = foOptD;
                bestI = i;
                bestJ = j;
            }
        }
    }

    double ofv = s.ofv + bestFoOpt;

    if (bestI != -1 && bestJ != -1 && ofv < s.ofv)
    {
        TSol sol = InsertNode(s, routeLength, bestI, bestJ);

        int deliveriesSize = deliveries.size();
        for (int i = 0; i < deliveriesSize; i++)
        {
            sol.vec[i + routeLength].sol = sortedDeliveriesByValue[i];
        }

        sol.ofv = ofv;
        sol.ofvR += bestFoOptR;
        sol.ofvD -= bestFoOptD;

        s = sol;
        // sanityCheck(s);
    }
}

void LS2Exact(TSol &s) {
    int bestI = -1, bestJ = -1;
    double bestFoOptD = 0, bestFoOptR = 0, bestFoOpt = 0;
    std::vector<int> bestDeliveries(deliveries.size());

    int routeLength = dist.size();

    for (int i = 1; i < routeLength; i++)
    {
        for (int j = (i + 1) % routeLength; j != (i + 1) % routeLength - 1 && j != i; j = (j + 1) % routeLength)
        {
            if (j == 0) {
                continue;
            }

            int vi  = s.vec[i].sol;
            int viP = s.vec[(i+1)%routeLength].sol;
            int viM = s.vec[i-1].sol;

            int vjM = s.vec[j - 1].sol;
            int vj  = s.vec[j].sol;
            int vjP = s.vec[(j+1)%routeLength].sol;

            double foOptR;

            if (i < j) {
                foOptR = - dist[vj][vjP]
                    - dist[viM][vi]
                    - dist[vi][viP]
                    + dist[vj][vi]
                    + dist[vi][vjP]
                    + dist[viM][viP];
            }
            else
            {
                foOptR = -dist[vjM][vj]
                    - dist[viM][vi]
                    - dist[vi][viP]
                    + dist[vi][vj]
                    + dist[vjM][vi]
                    + dist[viM][viP];

            }

            TSol newSol = InsertNode(s, routeLength, i, j);

            std::vector<int> optimumDeliveries = OptimumDeliveries(newSol);

            double foOptD = s.ofvD - SumDeliveries(optimumDeliveries);

            double foOpt = foOptR + foOptD;

            if (foOpt < bestFoOpt) {
                bestFoOpt = foOpt;
                bestFoOptR = foOptR;
                bestFoOptD = foOptD;
                bestI = i;
                bestJ = j;
                bestDeliveries = optimumDeliveries;
            }
        }
    }

    double ofv = s.ofv + bestFoOpt;

    if (bestI != -1 && bestJ != -1 && ofv < s.ofv)
    {
        TSol sol = InsertNode(s, routeLength, bestI, bestJ);

        int deliveriesSize = deliveries.size();
        for (int i = 0; i < deliveriesSize; i++)
        {
            sol.vec[i + routeLength].sol = bestDeliveries[i];
        }

        sol.ofv = ofv;
        sol.ofvR += bestFoOptR;
        sol.ofvD -= bestFoOptD;

        s = sol;

        // sanityCheck(s);
    }
}

void FreeMemoryProblem()
{
    //specific problem
    dist.clear();
    deliveries.clear();
}

void printSolution(TSol &s, int n)
{
    for (int i = 0; i < n; i++) {
        printf("%d, %d, %lf\n", i, s.vec[i].sol, s.vec[i].rk);
    }
}

void sanityCheck(TSol s) {
    int nodeSize = dist.size();
    int deliverySize = deliveries.size();

    double ofv = s.ofv;
    double ofvR = s.ofvR;
    double ofvD = s.ofvD;

    int copy[nodeSize + deliverySize];

    bool negative = false;
    for (int i = 0; i < nodeSize + deliverySize; i++) {
        copy[i] = s.vec[i].sol;

        if (s.vec[i].sol < 0) {
            negative = true;
        }
    }

    if (!negative) {
        SelectDeliveries(s, deliverySize + nodeSize, nodeSize);
    }
    CalculateFitness(s, nodeSize + deliverySize , nodeSize);

    if (s.ofv != ofv || s.ofvR != ofvR || s.ofvD != ofvD) {
        printf("\nERRO SANITY CHECK 1\n");
    }

    TSol before(s);
    auto optimumDeliveries = OptimumDeliveries(s);

    for (int i = 0; i < deliverySize; i++) {
        s.vec[i + nodeSize].sol = optimumDeliveries[i];
    }

    auto fitness = CalculateFitness(s, nodeSize + deliverySize, nodeSize);

    if (fitness > s.ofv) {
        printf("\nERRO SANITY CHECK 2\n");
        printf("Fitness: %f - Ofv %f\n", fitness, s.ofv);
        for (int i = 0; i < nodeSize + deliverySize; i++) {
            printf("%d ", s.vec[i].sol);
        }
        printf("\n");
        for (int i = 0; i < nodeSize + deliverySize; i++) {
            printf("%d ", before.vec[i].sol);
        }
        exit(1);
    }
}