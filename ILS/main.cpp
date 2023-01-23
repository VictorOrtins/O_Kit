#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <math.h>

/*
 typedef struct Solucao{
    vector<int> sequence;
    double valorObj;
 } Solucao;
Assim que é a struct de solução em C*/

/*
    CL = V\V' significa que CL é igual ao conjunto de vértices de V que não estão em V'
    CL, então, seria o complemento de V'
*/

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

        void pushBackSequence(int value){
            sequence.push_back(value);
        }

        void insertSequenceAt(int value, int position){

            if(position > sequence.size() - 1){
                return;
            }
            std::vector<int>::iterator it = sequence.begin() + position;
            sequence.insert(it, value);
        }

        void swapSequence(int position1, int position2){
            std::swap(sequence.at(position1), sequence.at(position2));
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
            printf("{ insertedNode: %d, removedEdge: %d, cost: %f }\n", insertedNode, removedEdgeIndex, cost);
        }
};

Solution iteratedLocalSearch(int maxIter, int maxIterILS);

Solution construction();

Solution perturbation(Solution bestSolution);

void localSearch(Solution* possibleSolution);

bool greaterCost(Solution first, Solution second);

void whileILS(Solution solution, int *iterILS);

std::vector<InsertionInfo> insertionCostCalculation(Solution &solution, std::vector<int> &complement);

std::vector<int> choose3RandomNodes();

std::vector<int> remainingNodes(std::vector<int> sequence);

void insertionSort(std::vector<InsertionInfo> &sequence);

void updateSolution(Solution &solution, InsertionInfo choosen, std::vector<int> &complement);

bool bestImprovementSwap(Solution *solution);


int main(void){

    for(int i = 0; i < numberOfVertices; i++){
        vertices.push_back(i + 1);
    }

    Solution solution = construction();

    printf("\n");
    solution.print();


    // Solution solution;
    // solution.pushBackSequence(vertices.at(0));
    // solution.pushBackSequence(vertices.at(1));
    // solution.pushBackSequence(vertices.at(2));

    // std::vector<int> complement;
    // for(int i = solution.getSequence().size(); i < numberOfVertices; i++){
    //     complement.push_back(i);
    // }

    // std::vector<InsertionInfo> insertionCost = insertionCostCalculation(solution, complement);

    // for(auto insertion: insertionCost){
    //     insertion.print();
    // }

    return 0;
}

Solution construction(){
    Solution solution;

    solution.setSequence(choose3RandomNodes());
    std::vector<int> complement = remainingNodes(solution.getSequence());

    double alpha;
    int choosen;

    while(!complement.empty()){
        std::vector<InsertionInfo> insertionCost = insertionCostCalculation(solution, complement);
        insertionSort(insertionCost);
        alpha = (double) rand() / RAND_MAX;
        choosen = rand() % ((int) ceil (alpha * insertionCost.size()));
        updateSolution(solution, insertionCost.at(choosen), complement);
    }

    return solution;
}

Solution perturbation(Solution bestSolution){
    Solution solution;
    return solution;
}

void localSearch(Solution* possibleSolution){

}

bool greaterCost(Solution first, Solution second){
    return first.getCost() < second.getCost();
}

Solution iteratedLocalSearch(int maxIter, int maxIterILS){
    Solution bestOfAll;
    int iterILS;

    bestOfAll.setCost(INFINITY);

    for(int i = 0; i < maxIter; i++){
        Solution solution = construction(); //Construção de uma solução baseado em "palpites educativos"
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
    int solutionSequenceSize = solution.getSequence().size();

    std::vector<InsertionInfo> insertionCost( (solutionSequenceSize) * (complement.size()));

    int l = 0;
    for(int a = 0, b = 1; a < solutionSequenceSize - 1; a++, b++){
        int vertexI = solution.getSequence().at(a);
        int vertexJ = solution.getSequence().at(b);

        for(auto vertexK : complement){
            insertionCost[l].setCost(edgeCost[vertexI][vertexK] + edgeCost[vertexJ][vertexK] - edgeCost[vertexI][vertexJ]);
            insertionCost[l].setInsertedNode(vertexK);
            insertionCost[l].setRemovedEdgeIndex(a);
            l++;
        }
    }

    return insertionCost;
}

std::vector<int> choose3RandomNodes(){
    std::vector<int> ret;

    for(int i = 0; i < 3; i++){
        ret.push_back( rand() % numberOfVertices);
    }

    return ret;
}

std::vector<int> remainingNodes(std::vector<int> sequence){
    std::vector<int> ret;

    for(int i = 0; i < numberOfVertices; i++){
        if(std::find(sequence.begin(), sequence.end(), i) == sequence.end()){
            ret.push_back(i);
        }
    }

    return ret;
}

void insertionSort(std::vector<InsertionInfo> &sequence){
    int j;
    InsertionInfo temp;
    for(int i = 1; i < sequence.size(); i++){
        temp = sequence.at(i);
        j = i - 1;
        while(j >= 0 && sequence.at(j).getCost() < temp.getCost()){
            sequence[j + 1] = sequence[j];
            j--;
        }
        sequence[j + 1] = temp;
    }
}

void updateSolution(Solution &solution, InsertionInfo choosen, std::vector<int> &complement){
    solution.insertSequenceAt(choosen.getInsertedNode(), choosen.getRemovedEdgeIndex() + 1);

    for(int i = 0; i < complement.size(); i++){
        if(complement.at(i) == choosen.getInsertedNode()){
            complement.erase(complement.begin() + i);
            break;
        }
    }
}

bool bestImprovementSwap(Solution *solution){
    double bestDelta = 0;
    double best_i, best_j;

    int i_value, i_value_next, i_value_prev;
    int j_value, j_value_next, j_value_prev;
    double delta;

    for(int i = 1; i < solution->getSequence().size() - 1; i++){
        i_value = solution->getSequence().at(i);
        i_value_next = solution->getSequence().at(i + 1);
        i_value_prev = solution->getSequence().at(i - 1);

        for(int j = i + 1; j < solution->getSequence().size() - 1; j++){
            j_value = solution->getSequence().at(j);
            j_value_next = solution->getSequence().at(j + 1);
            j_value_prev = solution->getSequence().at(j - 1);

            //Revisar essa fórmula do delta.
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
        solution->swapSequence(best_i, best_j);
        solution->setCost(solution->getCost() + delta);
        return true;
    }

    return false;
}





