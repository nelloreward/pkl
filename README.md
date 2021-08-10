# P(K)L sequence construction

This repo contains software to accompany the paper Arbitrary-length analogs to de Bruijn sequences, by [Abhinav Nellore](https://nellore.bio) and [Rachel Ward](https://sites.google.com/prod/view/rward).

# Quick start

* Construct a P(K)L sequence of length <sequence_length> on the alphabet {0, 1, ..., K - 1} of size <alphabet_size> with

        python3 pkl.py construct -l <sequence_length> -k <alphabet_size>
* Check that a sequence <sequence> is a P(K)L sequence with

        echo "<sequence>" | python3 pkl.py check
* Test that the software works at various sequence lengths and alphabet sizes with

        python3 pkl.py test
  
* Learn more about command-line options by running each of
  
        python3 pkl.py check --help
        python3 pkl.py construct --help
        python3 pkl.py test --help

 # Regenerating Table 1 from paper
  
 Run
  
         
