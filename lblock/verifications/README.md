# Experimental Verification

We have provided the required codes to verify the probability of our DL distinguishers experimentally. To specify the input/output differences, the number of rounds for DL distinguisher, as well as the number of DL queries open [`difflin.h`](difflin.h) and adjust the corresponding variables. Next, to compile the code and then perform the experimental evaluation, run the following commands:

```sh
make
./difflin 0
```