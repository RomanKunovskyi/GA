using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Term_Papper_GA
{
    public class Individuum : ICloneable
    {
        public Chromosome[][] Chromosomes { get; private set; }
        public double Eval { get; private set; }
        public Individuum(Chromosome[][] chromosomes, ValueContainer container)
        {
            Chromosomes = new Chromosome[chromosomes.Length][];
            for (int i = 0; i < chromosomes.Length; i++)
            {
                Chromosomes[i] = new Chromosome[chromosomes[i].Length];
                for (int j = 0; j < chromosomes[i].Length; j++)
                {
                    Chromosomes[i][j] = chromosomes[i][j];
                }
            }
            Eval = СalculateEval(container);
        }
        public Individuum(int a, int b, int n, int chromosomeLength, ValueContainer container)
        {
            Chromosomes = new Chromosome[n+1][];
            for (int i = 0; i < n + 1; i++)
            {
                Chromosomes[i] = new Chromosome[2*i+1];
                for (int j = 0; j < 2 * i + 1; j++)
                {
                    Chromosomes[i][j] = new Chromosome(chromosomeLength, a, b);
                }
            }
            Eval = СalculateEval(container);
        }

        private Individuum(Chromosome[][] chromosomes, double eval)
        {
            Chromosomes = new Chromosome[chromosomes.Length][];
            for (int i = 0; i < chromosomes.Length; i++)
            {
                Chromosomes[i] = new Chromosome[chromosomes[i].Length];
                for (int j = 0; j < chromosomes[i].Length; j++)
                {
                    Chromosomes[i][j] = (Chromosome)chromosomes[i][j].Clone();
                }
            }
            Eval = eval;
        }

        private double СalculateEval(ValueContainer container)
        {
            return container.GetEval(Chromosomes.Select(chromosome => chromosome.Select(arr => arr.Number).ToArray()).ToArray());
        }

        public override string ToString()
        {
            string str = "";
            int i = 0;
            foreach (var arr in Chromosomes)
            {
                int k = i;
                str += $"{k}: ";
                foreach (var chromosome in arr)
                {
                    str += $"{chromosome}, ";
                }
                str += "\n";
                i++;
            }
            return str + "Eval = " + Eval;
        }

        public object Clone()
        {
            return new Individuum(Chromosomes, Eval);
        }
    }
}
