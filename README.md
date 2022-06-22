## Construction of arbitrary-length analogs to de Bruijn sequences

This repo contains software for constructing $P^{(K)}_L$-sequences to accompany our [CPM 2022](http://www.stringology.org/event/CPM2022/) paper [*Arbitrary-length analogs to de Bruijn sequences*](https://drops.dagstuhl.de/opus/volltexte/2022/16136/), by [Abhinav Nellore](https://nellore.bio) and [Rachel Ward](https://sites.google.com/prod/view/rward).

## Citation

If you use our software, please cite this repo as well as our paper:

        @InProceedings{nellore_et_al:LIPIcs.CPM.2022.9,
          author = {Nellore, Abhinav and Ward, Rachel},
          title = {{Arbitrary-Length Analogs to de Bruijn Sequences}},
          booktitle = {33rd Annual Symposium on Combinatorial Pattern Matching (CPM 2022)},
          pages = {9:1--9:20},
          series = {Leibniz International Proceedings in Informatics (LIPIcs)},
          ISBN = {978-3-95977-234-1},
          ISSN = {1868-8969},
          year = {2022},
          volume = {223},
          editor = {Bannai, Hideo and Holub, Jan},
          publisher = {Schloss Dagstuhl -- Leibniz-Zentrum f{\"u}r Informatik},
          address = {Dagstuhl, Germany},
          URL =	{https://drops.dagstuhl.de/opus/volltexte/2022/16136},
          URN =	{urn:nbn:de:0030-drops-161361},
          doi =	{10.4230/LIPIcs.CPM.2022.9},
          annote = {Keywords: de Bruijn sequence, de Bruijn word, Lempel’s D-morphism, Lempel’s homomorphism}
        }

## Errata

The published version of our paper has a few minor errors resolved as follows:
* In the procedure `GeneratePKL` (Algorithm 2), to properly address de Bruijn sequence-length edge cases, use $N := \lfloor \log_K L \rfloor + 1$ instead of $N := \lceil \log_K L \rceil$, except in the call to `LiftAndJoin` (Algorithm 1), where $N := \lceil \log_K L \rceil$ should be used as in the text. References to $N$ in the proof of Theorem 9 should interpreted to respect this change.
* For the same reason, in the procedure `GenerateP2L` (Algorithm 3), use $N := \lfloor \log_2 L \rfloor + 1$ instead of $N := \lceil \log_2 L \rceil$.
* Equation (8) should read as $L_j = K \cdot L_{j-1} + d_j \quad j \in ${$1, 2, \ldots, N-1$}.

## Use at the command line

* Construct a $P^{(K)}_L$-sequence of length $L$ on the alphabet {0, 1, $\ldots$, $K$ - 1} with

        python3 pkl.py construct -l L -k K
  Running

        time python3 pkl.py construct -l 100000000 -k 7
  on a MacBook Air M1, 2020 with 16 GB of RAM using Python 3.7.3 (default, Jun 19 2019, 07:38:49) Clang 10.0.1 (clang-1001.0.46.4) on darwin gave

        real	1m9.755s
        user	0m57.393s
        sys	0m6.693s
* Check that a sequence <sequence> is a $P^{(K)}_L$-sequence with

        echo "<sequence>" | python3 pkl.py check
* Validate that the software successfully generates a $P^{(K)}_L$-sequence for any combination of 1 <= sequence length <= 10000 and 2 <= alphabet size <= 20 with

        python3 pkl.py test
  (We ran this command on a MacBook Air M1, 2020 using Python 3.7.3 (default, Jun 19 2019, 07:38:49) Clang 10.0.1 (clang-1001.0.46.4) on darwin and encountered no errors.)
        
* Learn more about command-line options by running each of
  
        python3 pkl.py check --help
        python3 pkl.py construct --help
        python3 pkl.py test --help

## Use in your software
        
        from pkl import pkl_via_lempels_lift
        
        pkl_sequence = pkl_via_lempels_lift(4, 1000) # Stores a length-1,000 $P^{(K)}_L$-sequence on the size-4 alphabet {0, 1, 2, 3} in the list pkl_sequence

## Regenerate Table 1 from paper
        
Run
        
        for i in {1..32}; do python3 count.py -k 2 -l $i; done
 
in Bash to perform exhaustive searches recovering the counts of binary $P^{(K)}_L$-sequences for sequence lengths between 1 and 32 displayed in Table 1 of the paper.
  
         
