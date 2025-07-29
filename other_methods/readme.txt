
1) x^2 case: we use the centred CI trick to use variance as a bound
- benchmark methods: Chebysev and Catoni
2) bounded: we compare to old methods, but most important are WSR, Emp Bernstein (Maurer and Pontil), Hoeffding
3) sub-gaussian: chernoff trick for CI

Folder called ~/cloud_results:
/bdd: has bdd_n.csv comparison where n = sample size
/x2: as x2_n.csv
/subg: as subg_n.csv