```sh
# Run the fit routine for DiRICH channels 960-991 and 1056-1087 with trigger at 994
./fit -i pulser_22265115533.A.root -d 960 1056 -t 994

# For fitting the histograms with gaus and plotting them
root pulser_22265115533.H.root

## number can be anything between 0-63
> .x fit.C(0)

# For plotting the histograms of sigmas, means etc
root pulser_22265115533.D.root
> .x plot.C
```

More info is there is `./fit -h`
