from __future__ import print_function, division
from argparse import ArgumentParser

def main(args):

    ## casa conversion task
    importuvfits(fitsfile=args.uvfits, vis="{}native_res.ms".format(args.obs))

    ## perform channel donw-sampling
    mstransform(vis="{}native_res.ms".format(args.obs),outputvis="{}.ms".format(args.obs),chanaverage=True, chanbin=(4))
    

if __name__ == "__main__":
    parser = ArgumentParser("uvfits2ms", description="script to convert uvfits to ms format using CASA tasks")
    parser.add_argument("--uvfits", required=True, help="the input uvfits file")
    parser.add_argument("--obs", required=True, help="the obsid of the observation")
    args = parser.parse_args()

    main(args)