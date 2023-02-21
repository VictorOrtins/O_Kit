#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <math.h>

std::vector<int> *solutionSequence;
const int numberOfVertices = 10;
std::vector<int> vertices;
int edgeCost[numberOfVertices][numberOfVertices];

class Solution{
    private:
        std::vector<int> sequence;
        double cost;

    public:
        Solution(std::vector<int> seq, double obj):sequence{seq}, cost{obj}{}

        Solution(){}

        void print(){
            for(int i = 0; i < sequence.size(); i++){
                std::cout << sequence.at(i);
                if(i != sequence.size() - 1){
                    std::cout << " -> ";
                }
            }
            std::cout << std::endl;
        }

        void setSequence(std::vector<int> seq){
            sequence = seq;
        }

        std::vector<int> getSequence(){
            return sequence;
        }

        void setCost(double obj){
            cost = obj;
        }

        double getCost(){
            return cost;
        }

        std::vector<int>* getSequencePointer(){
            return &sequence;
        }
};

class InsertionInfo{
    private:
        int insertedNode;
        int removedEdgeIndex;
        double cost;
    
    public:
        int getInsertedNode(){
            return insertedNode;
        }

        int getRemovedEdgeIndex(){
            return removedEdgeIndex;
        }

        double getCost(){
            return cost;
        }

        void setInsertedNode(int node){
            insertedNode = node;
        }

        void setRemovedEdgeIndex(int edge){
            removedEdgeIndex = edge;
        }

        void setCost(double value){
            cost = value;
        }

        void print(){
            printf("{ insertedNode: %d, removedEdgeIndex: %d, cost: %f }\n", insertedNode, removedEdgeIndex, cost);
        }
};


bool greaterCost(Solution first, Solution second);

std::vector<int> choose3RandomNodes();

std::vector<int> remainingNodes(std::vector<int> sequence);

void insertionSort(std::vector<InsertionInfo> &sequence);

bool vectorContains(std::vector<int> vec, int value);

bool bestImprovementSwap(Solution *solution);

bool bestImprovement2Opt(Solution *solution);

bool bestImprovementReinsertion(Solution *solution);

bool bestImprovementOrOpt(Solution *solution, int optNumber);

bool vectorContains(std::vector<int> vec, int value);

void orOptUpdate(std::vector<int> *solutionSequence, int opt, int best_i, int best_j);

Solution iteratedLocalSearch(int maxIter, int maxIterILS); 

Solution construction(); 

Solution perturbation(Solution bestSolution); 

void localSearch(Solution* possibleSolution); 

void whileILS(Solution solution, int *iterILS);

std::vector<InsertionInfo> insertionCostCalculation(Solution &solution, std::vector<int> &complement); 

void updateSolution(Solution &solution, InsertionInfo choosen, std::vector<int> &complement);

void printVector(std::vector<int> vec, int i_value);

void setEdgeCosts();

std::vector<InsertionInfo> insertionInfo();

double deltaUpdate(int opt, int i_value[], int j_value[]);

int getNextIndex(int currentIndex, int maxIndex);

int getSegmentSize(int nVertices);

bool validSegment(int segment_index, int segmentSize, std::vector<int> *solutionSequence);

bool swap_segments(std::vector<int> *vec, int start1, int size1, int start2, int size2);

int main(void){

    time_t t;
    // srand(time(&t));

    setEdgeCosts();

    for(int i = 0; i < numberOfVertices; i++){
        for(int j = 0; j < numberOfVertices; j++){
            printf("edgeCost[%d][%d] = %d\n", i, j, edgeCost[i][j]);
        }
        printf("\n");
    }

    for(int i = 0; i < numberOfVertices; i++){
        vertices.push_back(i + 1);
    }

    Solution solution = iteratedLocalSearch(50, 10);
    solution.print();
    printf("%f\n", solution.getCost());

    // bestImprovementSwap(&solution);

    // printf("\n");
    // solution.print();

    // bestImprovement2Opt(&solution);

    // printf("\n");
    // solution.print();

    // bestImprovementReinsertion(&solution);

    // printf("\n");
    // solution.print();

    // bestImprovementOrOpt(&solution, 2);

    // printf("\n");
    // solution.print();

    // bestImprovementOrOpt(&solution, 3);

    // printf("\n");
    // solution.print();

    // swap_segments(solution.getSequencePointer(), 1, 3, 5 , 9);

    // printf("\n");
    // printVector(solution.getSequence(), solution.getSequence().size());

    // solution = perturbation(solution);

    // printf("\n");
    // solution.print();


    return 0;
}

void setEdgeCosts(){
    for(int i = 0; i < numberOfVertices; i++){
        for(int j = i; j < numberOfVertices; j++){

            if(i == j){
                edgeCost[i][j] = 0;
                continue;
            }

            int num = rand() % 30 + 1;
            edgeCost[i][j] = num;
            edgeCost[j][i] = num;
        }
    }
}

std::vector<InsertionInfo> insertionInfo(){
    std::vector<InsertionInfo> vec;
    InsertionInfo temp;
    int num = rand() % 50 + 10;
    for(int i = 0; i < num; i++){
        temp.setCost( rand() % 30 - 10);
        temp.setInsertedNode( rand() % 10 + 1);
        temp.setRemovedEdgeIndex(rand() % 10 + 1);
        vec.push_back(temp);
    }


    return vec;
}

bool greaterCost(Solution first, Solution second){
    return first.getCost() < second.getCost();
}

std::vector<int> choose3RandomNodes(){
    std::vector<int> ret;

    ret.push_back(1);
    for(int i = 0; i < 3; i++){
        int vertex = (rand() % (numberOfVertices - 1)) + 2;
        if(vectorContains(ret, vertex) == true){
            i--;
            continue;
        }

        ret.push_back(vertex); //3 vértices aleatórios
    }

    ret.push_back(1);

    return ret;
}

std::vector<int> remainingNodes(std::vector<int> sequence){
    std::vector<int> ret;

    for(int i = 2; i < numberOfVertices + 1; i++){
        if(vectorContains(sequence, i) == true){
            continue;
        }

        ret.push_back(i);  //Apenas vértices que não estão
        //em sequence
    }

    return ret;
}

void insertionSort(std::vector<InsertionInfo> &sequence){
    int j;
    InsertionInfo temp;

    for(int i = 1; i < sequence.size(); i++){
        temp = sequence.at(i);
        j = i - 1;
        while(j >= 0 && sequence.at(j).getCost() > temp.getCost()){ //Ordem crescente
            sequence[j + 1] = sequence[j];
            j--;
        }
        sequence[j + 1] = temp;
    }
}

bool vectorContains(std::vector<int> vec, int value){

    if(std::find(vec.begin(), vec.end(), value) == vec.end()){ //não contém
        return false;
    }

    return true;
}

void printVector(std::vector<int> vec, int i_value){
    
    for(int i = 0; i < i_value; i++){
        std::cout << vec.at(i) << " ";
    }

    std::cout << std::endl;
}

bool bestImprovementSwap(Solution *solution){
    solutionSequence = solution->getSequencePointer();

    double bestDelta = 0; 
    int best_i, best_j;

    int i_value, i_value_next, i_value_prev;
    int j_value, j_value_next, j_value_prev;
    double delta;

    for(int i = 1; i < solutionSequence->size() - 1; i++){
        i_value = solutionSequence->at(i) - 1;
        i_value_next = solutionSequence->at(i + 1) - 1;
        i_value_prev = solutionSequence->at(i - 1) - 1;

        for(int j = i + 1; j < solutionSequence->size() - 1; j++){
            j_value = solutionSequence->at(j) - 1;
            j_value_next = solutionSequence->at(j + 1) - 1;
            j_value_prev = solutionSequence->at(j - 1) - 1 ;

            //Revisar essa fórmula do delta. REVISA ISSO AQUI!!
            delta = -edgeCost[i_value][i_value_prev] - edgeCost[i_value][i_value_next] + edgeCost[i_value_prev][i] + 
                    edgeCost[j_value][i_value_next] - edgeCost[j_value_prev][j_value] - edgeCost[j_value][j_value_next] + 
                    edgeCost[j_value][i_value] + edgeCost[i_value][j_value_next];

            if(delta < bestDelta){
                bestDelta = delta;
                best_i = i;
                best_j = j;
            } 
        }
    }

    if(bestDelta < 0){
        std::swap(solutionSequence->at(best_i), solutionSequence->at(best_j));
        solution->setCost(solution->getCost() + delta);
        return true;
    }

    return false;
}

bool bestImprovement2Opt(Solution *possibleSolution){
    solutionSequence = possibleSolution->getSequencePointer();

    double best_delta = 0, delta = 0;

    int best_i = 0, best_j = 0;
    int i_value = 0, i_next_value = 0;
    int j_value = 0, j_next_value = 0;

    for(int i = 0; i < solutionSequence->size() - 1; i++){
        i_value = solutionSequence->at(i) - 1;
        i_next_value = solutionSequence->at(i + 1) - 1;

        for(int j = i + 2; j < solutionSequence->size() - 1; j++){
            j_value = solutionSequence->at(j) - 1;
            j_next_value = solutionSequence->at(j + 1) - 1;

            delta = -edgeCost[i_value][i_next_value] - edgeCost[j_value][j_next_value] 
                    + edgeCost[i_value][j_value] + edgeCost[i_next_value][j_next_value];

            if(delta < best_delta){
                best_delta = delta;
                best_i = i;
                best_j = j;
            }
        }
    } 

    if(best_delta < 0){

        std::vector<int> subVector = {solutionSequence->begin() + best_i + 1, solutionSequence->begin() + best_j};

        solutionSequence->erase(solutionSequence->begin() + best_i + 1, solutionSequence->begin() + best_j);

        // //O(n)
        std::reverse(subVector.begin(), subVector.end());
        solutionSequence->insert(solutionSequence->begin() + best_i + 1, subVector.begin(), subVector.end());

        possibleSolution->setCost(possibleSolution->getCost() + delta);

        return true;
    }

    return false;
}

bool bestImprovementReinsertion(Solution *solution){
    return bestImprovementOrOpt(solution, 1); //Uma reinserção é só um Or-Opt de 1
}

bool bestImprovementOrOpt(Solution *possibleSolution, int opt){

    if(opt < 1 || opt > 3){
        return false;
    }

    solutionSequence = possibleSolution->getSequencePointer();
    const int solutionSequenceSize = solutionSequence->size();

    double best_delta = 0, delta = 0;

    int best_i = 0, best_j = 0;
    int i_value[5];
    // i_value_prev, i_value, i_value_next, i_value_2_next, i_value_3_next;

    int j_value[2];
    // int j_value, j_value_next;

    int j_next;
    int i_next[3];
    // int i_next, i_2_next, i_3_next;

    int j = 0;
    for(int i = 1; i < solutionSequenceSize - 1; i++){
        i_next[0] = getNextIndex(i, solutionSequenceSize - 1);

        i_value[1] = solutionSequence->at(i) - 1;
        i_value[2] = solutionSequence->at(i_next[0]) - 1;
        i_value[0] = solutionSequence->at(i - 1) - 1;

        if(i_value[1] == solutionSequence->at(0) - 1 || i_value[2] == solutionSequence->at(solutionSequenceSize - 1) - 1){
            continue;
        }

        j = i + 1;

        if(opt == 2 || opt == 3){
            i_next[1] = getNextIndex(i_next[0], solutionSequenceSize - 1);
            i_value[3] = solutionSequence->at(i_next[1]) - 1;
            j = getNextIndex(j, solutionSequenceSize - 1);

            if(i_value[3] == solutionSequence->at(0) - 1){
                continue;
            }

            if(opt == 3){
                i_next[2] = getNextIndex(i_next[1], solutionSequenceSize - 1);
                i_value[4] = solutionSequence->at(i_next[2]) - 1;
                j = getNextIndex(j, solutionSequenceSize - 1);
            } 
        }

        while(j != i - 1){

            if(j == 0 || j == solutionSequenceSize - 1){
                j = (j + 1) % solutionSequenceSize;                
                continue;
            }

            j_next = getNextIndex(j, solutionSequenceSize - 1);

            j_value[0] = solutionSequence->at(j) - 1;
            j_value[1] = solutionSequence->at(j_next) - 1;

            if(j_value[0] == solutionSequence->at(0) - 1 || j_value[0] == solutionSequence->at(solutionSequenceSize - 1) - 1){
                j = (j + 1) % solutionSequenceSize;                
                continue;                
            }

            delta = -edgeCost[i_value[0]][i_value[1]] - edgeCost[j_value[0]][j_value[1]] + edgeCost[i_value[1]][j_value[1]];

            delta += deltaUpdate(opt, i_value, j_value);

            if(delta < best_delta){
                best_delta = delta;
                best_i = i;
                best_j = j;
            }

            j = getNextIndex(j, solutionSequenceSize - 1);
        }
    }

    if(best_delta < 0){
        orOptUpdate(solutionSequence, opt, best_i, best_j);
        possibleSolution->setCost(possibleSolution->getCost() + best_delta);

        return true;
    }

    return false;
}

void orOptUpdate(std::vector<int> *solutionSequence, int opt, int best_i, int best_j){
    switch(opt){
        case 1:{
            int reinsertedValue = solutionSequence->at(best_i);
            solutionSequence->erase(solutionSequence->begin() + best_i);
            solutionSequence->insert(solutionSequence->begin() + best_j, reinsertedValue);
            break;
        }

        case 2:{
            std::vector<int> subVector {solutionSequence->begin() + best_i, solutionSequence->begin() + best_i + 2};
            solutionSequence->erase(solutionSequence->begin() + best_i, solutionSequence->begin() + best_i + 2);

            std::reverse(subVector.begin(), subVector.end());

            if(best_i < best_j){
                solutionSequence->insert(solutionSequence->begin() + best_j - 1, subVector.begin(), subVector.end());
            }
            else{
                solutionSequence->insert(solutionSequence->begin() + best_j + 1, subVector.begin(), subVector.end());
            }
            break;
        }

        case 3:{
            std::vector<int> subVector {solutionSequence->begin() + best_i, solutionSequence->begin() + best_i + 3};
            solutionSequence->erase(solutionSequence->begin() + best_i, solutionSequence->begin() + best_i + 3);

            std::reverse(subVector.begin(), subVector.end());

            if(best_i < best_j){
                solutionSequence->insert(solutionSequence->begin() + best_j - 3, subVector.begin(), subVector.end());
            }
            else{
                solutionSequence->insert(solutionSequence->begin() + best_j + 1, subVector.begin(), subVector.end());
            }
        }
    }
}

double deltaUpdate(int opt, int i_value[], int j_value[]){
    switch(opt){
        case 1:{
            return -edgeCost[i_value[1]][i_value[2]] + edgeCost[i_value[0]][i_value[2]] + edgeCost[j_value[0]][i_value[1]];
        }

        case 2:{
            return -edgeCost[i_value[2]][i_value[3]]
                    + edgeCost[i_value[0]][i_value[3]] + edgeCost[j_value[0]][i_value[2]];
        }

        case 3:{
            return -edgeCost[i_value[3]][i_value[4]] + edgeCost[j_value[0]][i_value[3]] + edgeCost[i_value[0]][i_value[4]];
        }
    }

    return 0.0;
}

Solution construction(){
    Solution solution;

    solution.setSequence(choose3RandomNodes()); //Escolha de 3 vértices aleatórios
    std::vector<int> complement = remainingNodes(solution.getSequence()); //Complemento dos 3 vértices

    double alpha;
    int choosen;

    while(!complement.empty()){ //Enquanto ainda tiverem vértices a colocar na solução

        std::vector<InsertionInfo> insertionCost = insertionCostCalculation(solution, complement); 
        //Todas possibilidades de inserção
        insertionSort(insertionCost); //Ordenar essas possibilidades em ordem crescente de custo

        alpha = (double) rand() / RAND_MAX; 
        choosen = rand() % ((int) ceil (alpha * insertionCost.size())); //Definir a inserção escolhida

        updateSolution(solution, insertionCost.at(choosen), complement); //Atualizar a solução, incluindo
        //remover um dos vértices do complemento
    }

    return solution;
}

Solution perturbation(Solution bestSolution){
    Solution solution = bestSolution;
    std::vector<int> *solutionSequence = bestSolution.getSequencePointer();
    int segmentSize;
    int segment1_init, segment2_init = 1, segment1_end = 0, segment2_end;

    while(segment1_end < segment2_init){
        segmentSize = getSegmentSize(solutionSequence->size());
        segment1_init = (rand() % (solutionSequence->size() - 2)) + 1;

        if(!validSegment(segment1_init, segmentSize, solutionSequence)){
            segment1_init = 1;
            segment1_end = 0;
            continue;
        }

        segment1_end = segment1_init + segmentSize;

        segment2_init = (rand() % (solutionSequence->size() - 2)) + 1;

        if(!validSegment(segment1_init, segmentSize, solutionSequence)){
            segment1_init = 1;
            segment1_end = 0;
            continue;
        }
    }

    swap_segments(solutionSequence, segment1_init, segment1_end - segment1_init + 1,
                 segment2_init, segment2_end - segment2_init + 1);

    return solution;
}

void localSearch(Solution* possibleSolution){
    std::vector<int> NL = {1,2,3,4,5}; //As 5 formas de melhorar a solução
    bool improved = false; //Se melhorou ou não
    while(NL.empty()){ //Enquanto ainda tiver forma de melhorar a solução
        int n = rand() % NL.size(); //A forma é escolhida aleatoriamente
        switch(NL[n]){
            case 1:
                improved = bestImprovementSwap(possibleSolution);
                break;

            case 2:
                improved = bestImprovement2Opt(possibleSolution);
                break;

            case 3:
                improved = bestImprovementReinsertion(possibleSolution);
                break;
            
            case 4:
                improved = bestImprovementOrOpt(possibleSolution, 2);
                break;

            case 5:
                improved = bestImprovementOrOpt(possibleSolution, 3);
                break;
        }

        if(improved){
            NL = {1,2,3,4,5}; //Se melhorou, tenta tudo de novo
        }
        else{
            NL.erase(NL.begin() + n); //Se não, continua com as opções ainda restantes
        }
    }
}

Solution iteratedLocalSearch(int maxIter, int maxIterILS){
    Solution bestOfAll;
    int iterILS;

    bestOfAll.setCost(INFINITY);

    for(int i = 0; i < maxIter; i++){
        Solution solution = construction(); //Construção de uma solução baseado em "palpites
        Solution best = solution;

        iterILS = 0;

        while(iterILS <= maxIterILS){ //Em tese, a local search em si fica aqui
            localSearch(&solution); //Tentando fazer melhoras na solução através de pequenas modificações

            if(greaterCost(solution, best)){
                best = solution;
                iterILS = 0; //Note que o iterILS vira zero quando uma solução melhor é encontrada
            }

            solution = perturbation(best); //Provoca uma "pertubação" na solução para melhorá-la
            iterILS++; 
        }

        if(greaterCost(best, bestOfAll)){
            bestOfAll = best;
        }
    }

    return bestOfAll;
}

std::vector<InsertionInfo> insertionCostCalculation(Solution &solution, std::vector<int> &complement){
    solutionSequence = solution.getSequencePointer();

    int solutionSequenceSize = solutionSequence->size();

    std::vector<InsertionInfo> insertionCost;

    InsertionInfo temp;
    //Calcula todas as possibilidades de inserção
    for(int a = 0, b = 1; a < solutionSequenceSize - 1; a++, b++){

        int vertexI = solutionSequence->at(a) - 1; //O vértice 7 tem as distâncias representadas na matrix na linha 6. Por isso o -1
        int vertexJ = solutionSequence->at(b) - 1;

        for(auto vertexK : complement){

            temp.setCost(edgeCost[vertexI][vertexK] + edgeCost[vertexJ][vertexK] - edgeCost[vertexI][vertexJ]);
            temp.setInsertedNode(vertexK);
            temp.setRemovedEdgeIndex(a);

            insertionCost.push_back(temp);
        }
    }

    return insertionCost;
}

void updateSolution(Solution &solution, InsertionInfo choosen, std::vector<int> &complement){
    solutionSequence = solution.getSequencePointer();

    solutionSequence->insert(solutionSequence->begin() + choosen.getRemovedEdgeIndex() + 1, choosen.getInsertedNode());
    //Colocar o vértice no lugar correspondente. se antes tinham os vértices i - j, e colocou-se um
    //k no meio (i - k - j), a função vai fazer isso ai

    for(int i = 0; i < complement.size(); i++){
        if(complement.at(i) == choosen.getInsertedNode()){ //Tirar o nó do complemento que foi colocado na solução
            complement.erase(complement.begin() + i);
            break;
        }
    }
}

int getNextIndex(int currentIndex, int maxIndex){
    return (currentIndex + 1) % (maxIndex);
}


int getSegmentSize(int nVertices){
    return 2 + ceil(rand() % ((nVertices)/10)); 
}

bool validSegment(int segment_index, int segmentSize, std::vector<int> *solutionSequence){
    return ((segment_index + segmentSize) % (solutionSequence->size() - 1)) > segment_index;
}

bool swap_segments(std::vector<int> *vec, int start1, int size1, int start2, int size2){

// Sort the segment indices to ensure that they are in increasing order
    if (start1 > start2) {
        std::swap(start1, start2);
        std::swap(size1, size2);
    }
    
    // Check that the segments do not overlap
    if (start1 + size1 > start2) {
        return false;
    }

    //Don't want to change segments that have the first or last element
    if(start1 == 0 || start2 == vec->size() - 1){
        return false;
    }

    if(start2 + size2 >= vec->size() - 1 || start1 + size1 >= vec->size() - 1){
        return false;
    }
    
    // Get iterators to the beginning and end of each segment
    auto seg1_begin = vec->begin() + start1;
    auto seg1_end = seg1_begin + size1;
    auto seg2_begin = vec->begin() + start2;
    auto seg2_end = seg2_begin + size2;
    
    // Swap the segments using std::rotate
    std::rotate(seg1_begin, seg2_begin, seg2_end);
    std::rotate(seg2_begin, seg1_end, seg1_end + size2);

    return true;
}







