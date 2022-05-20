using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Term_Papper_GA
{
    internal class GA
    {

        private readonly Random random = new Random();
        private readonly object locker = new object();

        private readonly int numberOfIndividuumInPopulation = 1000;
        private readonly int numberOfUselessPopulations = 50;
        private readonly int numberOfIndividuumInTournament = 2;

        private readonly double mutationProcent = 0.05;
        private readonly double crossoverProcent = 0.6;
        private readonly double chromosomeMutationProcent = 0.5;
        private readonly double chromosomeCrossoverProcent = 0.5;

        private readonly ValueContainer valueContainer;
        private readonly int chromosomeLength = 100;

        private readonly int n; //polinom length sum(1:n) sum(-m:m)
        private readonly int n_; //n'

        private readonly int A = -10;
        private readonly int B = 10;

        private readonly double E = 0.1;

        public GA()
        {
            n = 5;   //10 potolok   factorial problem
            n_ = 7;  //8 potolok   factorial problem 
            int number_x_theta = 5;
            int number_x_phi = 10;
            valueContainer = new ValueContainer(n, n_, number_x_theta, number_x_phi);

        }

        private List<Individuum> CreateFirstPopulation()
        {
            List<Individuum> population = new List<Individuum>(numberOfIndividuumInPopulation);
            Task[] tasks = new Task[numberOfIndividuumInPopulation];
            for (int i = 0; i < numberOfIndividuumInPopulation; i++)
            {
                tasks[i] = Task.Run(() =>
                    {
                        var individiuum = new Individuum(A, B, n, chromosomeLength, valueContainer);
                        lock (locker)
                        {
                            population.Add(individiuum);
                        }
                    });
            }
            Task.WaitAll(tasks);
            return population;
        }

        private Individuum GetTheBestIndividuum(List<Individuum> individuums)
        {
            Individuum bestIndividuum = individuums[0];
            for (int i = 1; i < individuums.Count; i++)
            {
                if (bestIndividuum.Eval > individuums[i].Eval)
                {
                    bestIndividuum = individuums[i];
                }
            }
            return (Individuum)bestIndividuum.Clone();
        }

        private Individuum Mutation(Individuum individuum)
        {
            Chromosome[][] chromosomes = new Chromosome[individuum.Chromosomes.Length][];
            for (int i = 0; i < chromosomes.Length; i++)
            {
                chromosomes[i] = new Chromosome[individuum.Chromosomes[i].Length];
                for (int j = 0; j < chromosomes[i].Length; j++)
                {

                    if (random.NextDouble() < chromosomeMutationProcent)
                    {
                        chromosomes[i][j] = new Chromosome(Mutation(individuum.Chromosomes[i][j].Array), A, B);
                    }
                    else
                    {
                        chromosomes[i][j] = (Chromosome)individuum.Chromosomes[i][j].Clone();
                    }

                }
            }
            return new Individuum(chromosomes, valueContainer);
        }
        private BitArray Mutation(BitArray arr)
        {
            BitArray newArr = (BitArray)arr.Clone();
            int randomNumber = random.Next() % arr.Length;
            newArr[randomNumber] = !newArr[randomNumber];
            return newArr;
        }

        //Crossover
        private (Individuum individuum1, Individuum individuum2) Crossover(Individuum individuum1, Individuum individuum2)
        {
            Chromosome[][] chromosomes1 = new Chromosome[individuum1.Chromosomes.Length][];
            Chromosome[][] chromosomes2 = new Chromosome[individuum2.Chromosomes.Length][];

            for (int i = 0; i < chromosomes1.Length; i++)
            {
                chromosomes1[i] = new Chromosome[individuum1.Chromosomes[i].Length];
                chromosomes2[i] = new Chromosome[individuum2.Chromosomes[i].Length];
                for (int j = 0; j < chromosomes1[i].Length; j++)
                {

                    if (random.NextDouble() < chromosomeCrossoverProcent)
                    {
                        (BitArray arr1, BitArray arr2) = Crossover(
                        individuum1.Chromosomes[i][j].Array, individuum2.Chromosomes[i][j].Array);

                        chromosomes1[i][j] = new Chromosome(arr1, A, B);
                        chromosomes2[i][j] = new Chromosome(arr2, A, B);
                    }
                    else
                    {
                        chromosomes1[i][j] = (Chromosome)individuum1.Chromosomes[i][j].Clone();
                        chromosomes2[i][j] = (Chromosome)individuum2.Chromosomes[i][j].Clone();
                    }
                }
            }
            return (new Individuum(chromosomes1, valueContainer), new Individuum(chromosomes2, valueContainer));
        }
        private (BitArray arr1, BitArray arr2) Crossover(BitArray arr1, BitArray arr2)
        {
            BitArray newArr1 = (BitArray)arr1.Clone();
            BitArray newArr2 = (BitArray)arr2.Clone();
            int randomNumber = random.Next(arr1.Length);
            for (int i = 0; i < randomNumber + 1; i++)
            {
                newArr1[i] = arr2[i];
                newArr2[i] = arr1[i];
            }
            return (newArr1, newArr2);
        }

        private void ConsolePrint(int populationNumber, Individuum bestIndividuumInPopulation, Individuum bestIndividuum)
        {
            string str = populationNumber + "_______________________________________\n" +
            "Best Individuum In Population: " +
            bestIndividuumInPopulation.Eval + "\n" +
            "Best Individuum: " +
            bestIndividuum.Eval +
            "\n" + "\n";

            Console.WriteLine(str);
        }
        private void ConsolePrintResult(int bestIndividuumPopulationNumber, Individuum bestIndividuum)
        {
            string str = bestIndividuumPopulationNumber + "____________________Result_____________________\n" +
            "Best Individuum: \n" +
            bestIndividuum.ToString() +
            "\n";

            Console.WriteLine(str);
        }
        
        public Individuum Run()
        {
            string forPlotBestInPopulation = "";
            string forPlotBest = "";

            int populationNumber = 0;
            int bestIndividuumPopulationNumber = populationNumber;
            List<Individuum> population = CreateFirstPopulation();
            Individuum bestIndividiuum = GetTheBestIndividuum(population);
            Individuum bestIndividiuumInPopulation = bestIndividiuum;
            Random random = new Random();

            Task[] tasks = new Task[numberOfIndividuumInPopulation];

            while (numberOfUselessPopulations >= populationNumber - bestIndividuumPopulationNumber || E < bestIndividiuum.Eval)
            {

                forPlotBestInPopulation += bestIndividiuumInPopulation.Eval + " ";
                forPlotBest += bestIndividiuum.Eval + " ";
                ConsolePrint(populationNumber, bestIndividiuumInPopulation, bestIndividiuum);

                List<Individuum> newPopulation = new List<Individuum>(numberOfIndividuumInPopulation);

                for (int i = 0; i < numberOfIndividuumInPopulation; i++)
                {
                    List<Individuum> tournamentParticipants = new List<Individuum>(numberOfIndividuumInTournament);
                    for (int j = 0; j < numberOfIndividuumInTournament; j++)
                    {
                        tournamentParticipants.Add(population[random.Next(population.Count)]);
                    }
                    newPopulation.Add(GetTheBestIndividuum(tournamentParticipants));
                }

                int partnerIndex = -1;
                for (int i = 0; i < numberOfIndividuumInPopulation; i++)
                {
                    int index = i;
                    tasks[i] = Task.Run(() =>
                    {
                        if (mutationProcent >= random.NextDouble())
                        {
                            newPopulation[index] = Mutation(newPopulation[index]);
                        }
                        if (crossoverProcent >= random.NextDouble())
                        {
                            int pIndex;
                            lock (locker)
                            {
                                if (partnerIndex == -1)
                                {
                                    partnerIndex = index;
                                    return;
                                }
                                pIndex = partnerIndex;
                                partnerIndex = -1;
                            }
                            (newPopulation[pIndex], newPopulation[index]) = Crossover(newPopulation[pIndex], newPopulation[index]);
                        }
                    });
                }
                Task.WaitAll(tasks);


                population = newPopulation;
                populationNumber += 1;
                bestIndividiuumInPopulation = GetTheBestIndividuum(population);
                if (bestIndividiuum.Eval > bestIndividiuumInPopulation.Eval)
                {
                    bestIndividiuum = bestIndividiuumInPopulation;
                    bestIndividuumPopulationNumber = populationNumber;
                }
            }

            using (StreamWriter sw = new StreamWriter("forPlotBestInPopulation.txt", false, Encoding.Default))
            {
                sw.Write(forPlotBestInPopulation.Replace(",", "."));
            }
            using (StreamWriter sw = new StreamWriter("forPlotBest.txt", false, Encoding.Default))
            {
                sw.Write(forPlotBest.Replace(",", "."));
            }

            ConsolePrintResult(bestIndividuumPopulationNumber, bestIndividiuum);
            Console.WriteLine("___________________________________________");        


            var error = valueContainer.Test(n_, bestIndividiuum.Chromosomes.Select(c => c.Select(c => c.Number).ToArray()).ToArray());

            foreach (var c in error)
            {
                Console.WriteLine(c);
            }

            return bestIndividiuum;
        }
    }
}
