//
//  main.cpp
//  4470GA
//
//  Worked by Group 7
//

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>
#include <math.h>

#define CLOCKS_PER_MS (clock_t(1000))

using namespace std;

//const string GENES = "0123456789";
//const string GENES = "0123456789abcdefghijklmnopqrstuvwxyz";
const string GENES = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ , .-;:_!\"#%&/()=?@${[]}";
const string TARGET = "3590541721";
//const int TARGET_LENGTH = 10;

struct Individual {
    string chromosome;
    double fitness;
};

double fRand(double fMin, double fMax)
{
    double f = static_cast<double>(rand()) / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

char randomGenes() {
    int num = rand() % GENES.size();
    return GENES[num];
}

string createGnome() {
    string gnome = "";
    for(int i = 0; i < TARGET.length(); i++)
        gnome += randomGenes();
    return gnome;
}

double fitness(string chromosome) {
    double value = 0.1;
    for(int i = 0; i < TARGET.length(); i++)
    {
        if(chromosome[i] == TARGET[i])
            value *= 10;
    }
    return value;
}

string singlePoint() {
    string mask = "";
    for(int i = 0; i < TARGET.length(); i++) {
        if (i < TARGET.length() / 2)
            mask += '1';
        else
            mask += '0';
    }
    return mask;
}

void crossover(Individual p1, Individual p2, vector<Individual>& newGen) {
    string mask = singlePoint();
    string newChromosome1 = "";
    string newChromosome2 = "";
    for(int i = 0; i < TARGET.length(); i++) {
        if (mask[i] == '1') {
            newChromosome1 += p1.chromosome[i];
            newChromosome2 += p2.chromosome[i];
        }   else {
            newChromosome2 += p1.chromosome[i];
            newChromosome1 += p2.chromosome[i];
        }
    }
    Individual newOffspring1 {newChromosome1, 0.0}; // fitness value will be calculated in step 5
    Individual newOffspring2 {newChromosome2, 0.0};
    newGen.push_back(newOffspring1);
    newGen.push_back(newOffspring2);
}

void mutateGenes(Individual& genes) {
    int mutatePosition = rand() % genes.chromosome.length();
    int num = rand() % GENES.size();
    genes.chromosome[mutatePosition] = GENES[num];
}

int main(int argc, const char * argv[]) {

    clock_t start;
    clock_t end;
    long solveTime;  //in milliseconds
    srand(static_cast<unsigned>(time(NULL)));  //  using the time seed from srand explanation
    
    const int POPULATION_SIZE = 100; //the number of hypotheses to be include in the population
    int generation = 0; //record the n'th generation
    double maxFitness = 0.0;
    double totalFitness = 0.0;
    double crossoverRate = 0.4;  //the fraction of the propulation to replace by Crossover
    double mutationRate = 0.1;  //the mutation rate
    int childNum = crossoverRate * POPULATION_SIZE; //the number of new offspring
    int parentNum = POPULATION_SIZE - childNum;  //the number survived in the old generation
    vector<Individual> population;
    double threshold = 0.999 * pow(10.0, static_cast<double>(TARGET.size() - 1));
    start = clock();  //start clocking
    // create initial population
    for(int i = 0; i < POPULATION_SIZE; i++)
    {
        string gnome = createGnome();
        double fitnessValue = fitness(gnome);  //compute Fitness(h)
        totalFitness += fitnessValue;
        Individual newIndividual {gnome, fitnessValue};
        population.push_back(newIndividual);
        if (fitnessValue > maxFitness) {
            maxFitness = fitnessValue;
        }
    }
    // loop through every new generation, set the Target length as the threshold
    while (maxFitness < threshold) {
        generation++;
        string bestGene;
        vector<Individual> newPopulation;
        int j;
        int k = 0;
        // step 1: probabilitiscally select members from the population
        while (k < population.size() && newPopulation.size() < parentNum) {
            if (k == (population.size() - parentNum + newPopulation.size()) || population[k].fitness / totalFitness > fRand(0, 1)) {
                newPopulation.push_back(population[k]);
            }
            k++;
        }
        vector<Individual> parents;
        // step 2: Crossover
        k = 0;
        // probabilitiscally select parents from the population
        while (k < population.size() && parents.size() < childNum) {
            if (k == (population.size() - childNum + parents.size()) || population[k].fitness / totalFitness > fRand(0, 1)) {
                parents.push_back(population[k]);
            }
            k++;
        }
        // reproduce offsprings and put them into the new population
        j = childNum - 1;
        k = 0;
        while (k < j) {
            crossover(parents[k], parents[j], newPopulation);
            k++;
            j--;
        }
        // step 3: Mutate
        set<int> mutateIndex;
        while (mutateIndex.size() < mutationRate * newPopulation.size()) {
            mutateIndex.insert(rand() % newPopulation.size());
        }
        for (int index : mutateIndex) {
            mutateGenes(newPopulation[index]);
        }
        // step 4: update
        population = newPopulation;
        // step 5: evalute the new population
        maxFitness = 0.0;
        totalFitness = 0.0;
        for(int i = 0; i < population.size(); i++)
        {
            population[i].fitness = fitness(population[i].chromosome);  //compute Fitness(h)
            totalFitness += population[i].fitness;
            if (population[i].fitness > maxFitness) {
                maxFitness = population[i].fitness;
                bestGene = population[i].chromosome;
            }
        }
        cout << "Generation: " << generation << "  " << bestGene << "  fitness: " << maxFitness << endl;
    }
    end = clock();
    solveTime = (end - start) / CLOCKS_PER_MS;
    cout << "\nTotal time (in millisecond): " << solveTime << endl;
    return 0;
}
