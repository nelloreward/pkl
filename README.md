# P(K)L-sequence construction

This repo contains software to accompany the paper Arbitrary-length analogs to de Bruijn sequences, by [Abhinav Nellore](https://nellore.bio) and [Rachel Ward](https://sites.google.com/prod/view/rward).

## Use at the command line

* Construct a P(K)L-sequence of length <sequence_length> on the alphabet {0, 1, ..., K - 1} of size <alphabet_size> with

        python3 pkl.py construct -l <sequence_length> -k <alphabet_size>
  Running

        time python3 pkl.py construct -l 100000000 -k 7
  on a MacBook Air M1, 2020 with 16 GB of RAM using Python 3.7.3 (default, Jun 19 2019, 07:38:49) Clang 10.0.1 (clang-1001.0.46.4) on darwin gave

        real    1m12.158s
        user    0m56.152s
        sys     0m8.044s
* Check that a sequence <sequence> is a P(K)L-sequence with

        echo "<sequence>" | python3 pkl.py check
* Validate that the software successfully generates a P(K)L-sequence for any combination of 1 <= sequence length <= 10000 and 2 <= alphabet size <= 20 with

        python3 pkl.py test
  (We ran this command on a MacBook Air M1, 2020 using Python 3.7.3 (default, Jun 19 2019, 07:38:49) Clang 10.0.1 (clang-1001.0.46.4) on darwin and encountered no errors.)
        
* Learn more about command-line options by running each of
  
        python3 pkl.py check --help
        python3 pkl.py construct --help
        python3 pkl.py test --help

## Use in your software
        
        from pkl import pkl_via_lempels_lift
        
        pkl_sequence = pkl_via_lempels_lift(4, 1000) # Stores a length-1,000 P(K)L sequence on the size-4 alphabet {0, 1, 2, 3} in the list pkl_sequence

## Regenerate Table 1 from paper
        
Run
        
        for i in {1..32}; do python3 count.py -k 2 -l $i; done
 
in Bash to perform exhaustive searches recovering the counts of binary P(K)L sequences for sequence lengths between 1 and 32 displayed in Table 1 of the paper.
  
         
