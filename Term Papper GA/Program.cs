using LegendreFunction;
using System.Collections;
using System.Diagnostics;
using Term_Papper_GA;
using MathNet.Numerics.RootFinding;
using MathNet.Numerics;


////test Legandre  legendre(N, X) in octaive
//var leg1 = Legendre.ALegendFnmSeri(3, -0.1);
//var leg2 = Legendre.ALegendFnmSeri(6, 0.3);
//var leg3 = Legendre.ALegendFnmSeri(8, -0.7);
//foreach (var l in leg1)
//{10
//    Console.WriteLine(l);
//}
//Console.WriteLine();
//foreach (var l in leg2)
//{
//    Console.WriteLine(l);
//}
//Console.WriteLine();
//foreach (var l in leg3)
//{
//    Console.WriteLine(l);
//}

GA ga = new GA();
ga.Run();

