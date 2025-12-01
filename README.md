# Nine and Ten Lonely Runners Code

This repository accompanies the article "[Nine and Ten Lonely Runners](https://arxiv.org/abs/2511.22427)" by Tanupat (Paul) Trakulthongchai. Our code is developed based on the [code](https://gite.lirmm.fr/mrosenfeld/the-lonely-runner-conjecture) provided by Matthieu Rosenfeld.

## Files
``lrc_for_nine_runners.cpp`` is the main code that checks whether $I(k,k+1,p)$ is empty for $k=8$.

``meta_lrc_nine.sh`` is the bash script that runs ``lrc_for_nine_runners.cpp`` for a fixed $k$ and a list of primes.

``results_nine.txt`` is the receipt of our own run that allows us to conclude in Section 5 of the paper that $I(k,k+1,p)$ is empty for all $p\in S_k$ for $k=8$. 

The same applies for 10 runners, if we replace ``nine`` with ``ten`` in the file names (and $k=8$ with $k=9$).


## Verifying our results
To verify our results for 9 runners, run the following bash script:

```bash
chmod +x meta_lrc_nine.sh
./meta_lrc_nine.sh lrc_for_nine_runners.cpp 8 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173 179 181 191 193 197 199 211 223 227 229 233 239 241 results_nine.txt
```

For 10 runners, run:
```bash
chmod +x meta_lrc_ten.sh
./meta_lrc_ten.sh lrc_for_ten_runners.cpp 9 137 139 149 151 157 163 167 173 179 181 191 193 197 199 211 223 227 229 233 239 241 251 257 263 269 271 277 281 283 293 307 311 313 317 331 337 347 349 353 359 367 373 379 383 389 397 401 results_ten.txt
```

## License

This code is released under the [MIT](https://choosealicense.com/licenses/mit/) license.
