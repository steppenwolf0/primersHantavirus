using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace WindowsFormsApp1
{
    
    public class PrimerUtils
    {
        // ✅ Complement mapping
        static readonly Dictionary<char, char> compMap = new Dictionary<char, char>
        {
            {'A','T'}, {'T','A'}, {'G','C'}, {'C','G'},
            {'N','N'}   // ✅ add this
        };

        static string ReverseComplement(string seq)
        {
            char[] rc = seq
                .Reverse()
                .Select(b => compMap.ContainsKey(b) ? compMap[b] : 'N') // ✅ safe
                .ToArray();

            return new string(rc);
        }


        // ✅ Tm computation
        public static double ComputeTm(string seq,
            double dna_nM = 50.0,
            double mv = 50.0,
            double dv = 1.5,
            double dntp = 0.6,
            double dmso = 0.0,
            double dmso_fact = 0.6,
            double formamide = 0.0)
        {
            seq = seq.ToUpper();
            int N = seq.Length;

            if (N < 2)
                throw new Exception("Sequence too short");

            if (dv < dntp)
                dv = dntp;

            double mv_eff = mv + 120.0 * Math.Sqrt(dv - dntp);

            var DH = new Dictionary<string, int>
        {
            {"AA",79},{"TT",79},{"AT",72},{"TA",72},
            {"CA",85},{"TG",85},{"GT",84},{"AC",84},
            {"CT",78},{"AG",78},{"GA",82},{"TC",82},
            {"CG",106},{"GC",98},{"GG",80},{"CC",80}
        };

            var DS = new Dictionary<string, int>
        {
            {"AA",222},{"TT",222},{"AT",204},{"TA",213},
            {"CA",227},{"TG",227},{"GT",224},{"AC",224},
            {"CT",210},{"AG",210},{"GA",222},{"TC",222},
            {"CG",272},{"GC",244},{"GG",199},{"CC",199}
        };

            double dh = 0;
            double ds = 0;

            // symmetry
            string rc_seq = ReverseComplement(seq);
            if (seq == rc_seq)
                ds += 14;

            // terminal corrections
            if ("AT".Contains(seq[0]))
            {
                ds -= 41;
                dh -= 23;
            }
            else
            {
                ds += 28;
                dh -= 1;
            }

            if ("AT".Contains(seq[N - 1]))
            {
                ds -= 41;
                dh -= 23;
            }
            else
            {
                ds += 28;
                dh -= 1;
            }

            // nearest-neighbor
            for (int i = 0; i < N - 1; i++)
            {
                string pair = seq.Substring(i, 2);

                if (DH.ContainsKey(pair) && DS.ContainsKey(pair)) // ✅ FIX
                {
                    dh += DH[pair];
                    ds += DS[pair];
                }
            }

            double delta_H = dh * -100.0;
            double delta_S = ds * -0.1;

            delta_S += 0.368 * (N - 1) * Math.Log(mv_eff / 1000.0);

            double R = 1.987;

            double tm_K = delta_H / (delta_S + R * Math.Log(dna_nM / 4e9));
            double tm_C = tm_K - 273.15;

            tm_C -= dmso * dmso_fact;

            double gc = GCContent(seq);
            tm_C += (0.453 * gc - 2.88) * formamide;

            return tm_C;
        }

        public static double GCContent(string seq)
        {
            int count = seq.Count(b => b == 'G' || b == 'C');
            return (double)count / seq.Length;
        }

        public static int GCClamp(string seq, int window = 5)
        {
            return seq.Substring(seq.Length - window)
                      .Count(b => b == 'G' || b == 'C');
        }

        public static double SelfAny(string seq)
        {
            string rc = ReverseComplement(seq);
            int n = seq.Length;
            int best = 0;

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    int k = 0;
                    while (i + k < n && j + k < n &&
                           seq[i + k] == rc[j + k])
                        k++;

                    if (k > best)
                        best = k;
                }
            }

            int score = best * 100;
            return score < 300 ? 0.0 : score / 100.0;
        }

        public static double SelfEnd(string seq, int k = 5)
        {
            string sub = seq.Substring(seq.Length - k);
            string rc = ReverseComplement(sub);

            int best = 0;

            for (int i = 0; i < k; i++)
            {
                for (int j = 0; j < k; j++)
                {
                    int m = 0;
                    while (i + m < k && j + m < k &&
                           sub[i + m] == rc[j + m])
                        m++;

                    if (m > best)
                        best = m;
                }
            }

            int score = best * 100;
            return score < 200 ? 0.0 : score / 100.0;
        }

        public static int Hairpin(string seq, int minLoop = 3, int minStem = 4)
        {
            int n = seq.Length;
            int best = 0;

            for (int i = 0; i < n; i++)
            {
                for (int j = i + minLoop + minStem; j < n; j++)
                {
                    int k = 0;
                    while (i - k >= 0 && j + k < n &&
                           seq[i - k] == compMap[seq[j + k]])
                    {
                        k++;
                    }

                    if (k >= minStem)
                        best = Math.Max(best, k);
                }
            }

            return best;
        }

        public static double EndStability(string seq)
        {
            var DG = new Dictionary<string, double>
        {
            {"AA",-1.0},{"TT",-1.0},{"AT",-0.88},{"TA",-0.58},
            {"CA",-1.45},{"TG",-1.45},{"GT",-1.44},{"AC",-1.44},
            {"CT",-1.28},{"AG",-1.28},{"GA",-1.30},{"TC",-1.30},
            {"CG",-2.17},{"GC",-2.24},{"GG",-1.84},{"CC",-1.84}
        };

            string tail = seq.Substring(seq.Length - 5);
            double dg = 0;

            for (int i = 0; i < tail.Length - 1; i++)
            {
                string pair = tail.Substring(i, 2);
                if (DG.ContainsKey(pair))
                    dg += DG[pair];
            }

            return Math.Abs(dg) * 0.6;
        }

        public static double Primer3Penalty(
            int length,
            double tm,
            double? gc = null,
            double posPenalty = 0.0,
            int sizeOpt = 20,
            double tmOpt = 60.0,
            double gcOpt = 50.0,
            double wSizeLt = 1.0,
            double wSizeGt = 1.0,
            //double wTmLt = 1.0, Original
            double wTmLt = 1.0,

            double wTmGt = 1.0,
            double wGcLt = 0.0,
            double wGcGt = 0.0,
            double wPos = 1.0)
        {
            double Linear(double value, double opt, double wLt, double wGt)
            {
                if (value < opt)
                    return (opt - value) * wLt;
                else
                    return (value - opt) * wGt;
            }

            double total = 0.0;

            total += Linear(length, sizeOpt, wSizeLt, wSizeGt);
            total += Linear(tm, tmOpt, wTmLt, wTmGt);

            if (gc.HasValue)
                total += Linear(gc.Value, gcOpt, wGcLt, wGcGt);

            total += wPos * posPenalty;

            return total;
        }
    }
}
