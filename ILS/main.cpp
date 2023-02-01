#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <math.h>
#include <bits/stdc++.h>
#include <stdio.h>

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
std::vector<int> *solutionSequence;

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

bool bestImprovement2Opt(Solution *solution);

bool bestImprovementReinsertion(Solution *solution);

bool bestImprovementOrOpt(Solution *solution, int optNumber);


int main(void){

    for(int i = 0; i < numberOfVertices; i++){
        vertices.push_back(i + 1);
    }

    Solution solution = construction();

    printf("\n");
    solution.print();

    bestImprovementSwap(&solution);

    printf("\n");
    solution.print();

    // bestImprovement2Opt(&solution);

    // printf("\n");
    // solution.print();



    // Solution solution;
    // std::vector<int> *sequence = solution.getSequencePointer();
    // sequence->push_back(vertices.at(0));
    // sequence->push_back(vertices.at(1));
    // sequence->push_back(vertices.at(2));

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
    Solution solution;
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
    solutionSequence = solution.getSequencePointer();

    int solutionSequenceSize = solutionSequence->size();

    std::vector<InsertionInfo> insertionCost( (solutionSequenceSize) * (complement.size()));

    int l = 0;
    //Calcula todas as possibilidades de inserção
    for(int a = 0, b = 1; a < solutionSequenceSize - 1; a++, b++){
        int vertexI = solutionSequence->at(a);
        int vertexJ = solutionSequence->at(b);

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
        ret.push_back( rand() % numberOfVertices); //3 vértices aleatórios
    }

    return ret;
}

std::vector<int> remainingNodes(std::vector<int> sequence){
    std::vector<int> ret;

    for(int i = 0; i < numberOfVertices; i++){
        if(std::find(sequence.begin(), sequence.end(), i) == sequence.end()){ //Apenas vértices que não estão
        //em sequências
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

bool bestImprovementSwap(Solution *solution){
    solutionSequence = solution->getSequencePointer();

    double bestDelta = 0; 
    int best_i, best_j;

    int i_value, i_value_next, i_value_prev;
    int j_value, j_value_next, j_value_prev;
    double delta;

    for(int i = 1; i < solutionSequence->size() - 1; i++){
        i_value = solutionSequence->at(i);
        i_value_next = solutionSequence->at(i + 1);
        i_value_prev = solutionSequence->at(i - 1);

        for(int j = i + 1; j < solutionSequence->size() - 1; j++){
            j_value = solutionSequence->at(j);
            j_value_next = solutionSequence->at(j + 1);
            j_value_prev = solutionSequence->at(j - 1);

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

bool bestImprovementReinsertion(Solution *solution){
    return bestImprovementOrOpt(solution, 1); //Uma reinserção é só um Or-Opt de 1
}

bool bestImprovement2Opt(Solution *possibleSolution){
    solutionSequence = possibleSolution->getSequencePointer();

    double best_delta = 0, delta = 0;

    int best_i = 0, best_j = 0;
    int i_value = 0, i_next_value = 0;
    int j_value = 0, j_next_value = 0;

    for(int i = 0; i < solutionSequence->size() - 1; i++){
        i_value = solutionSequence->at(i);
        i_next_value = solutionSequence->at(i + 1);

        for(int j = i + 2; j < solutionSequence->size() - 1; j++){
            j_value = solutionSequence->at(j);
            j_next_value = solutionSequence->at(j + 1);

            delta = -edgeCost[i][i + 1] - edgeCost[j][j + 1] + edgeCost[i][j] + edgeCost[i + 1][j + 1];

            if(delta < best_delta){
                best_delta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }

    if(best_delta < 0){
        std::vector<int> subVector = {solutionSequence->begin() + best_i + 1, solutionSequence->begin() + best_j};

        solutionSequence->erase(solutionSequence->begin() + best_i + 1, solutionSequence->end() + best_j);

        //O(n)
        std::reverse(subVector.begin(), subVector.end());
        solutionSequence->insert(solutionSequence->begin() + best_i + 1, subVector.begin(), subVector.end());

        return true;
    }

    return false;
}

bool bestImprovementOrOpt(Solution *possibleSolution, int opt){
    return false;
}