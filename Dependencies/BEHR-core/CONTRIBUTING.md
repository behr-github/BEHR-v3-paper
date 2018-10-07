# Contributing

If you want to contribute code to BEHR, we recommend that you open an issue on the relevant
repository labeled `proposed enhancement` and wait for feedback from one of the developers
before sinking significant time into the coding. You should explicitly discuss coauthorship
in future BEHR papers and be sure you are satisfied with the agreement before proceeding.

## Getting started

To submit your code, follow these steps for this or any of the dependencies:

1. Be sure you have a GitHub account
1. Fork the repository to your GitHub
1. Clone your fork to the computer you will do the work on.
1. Create a new branch off of the `develop` branch named `feature/<name>`, or
   `fix/<name>`, e.g. `fix/inconsistent-fill-values`. 
1. Once the new branch is done, pull `develop` from the BEHR repository and
   merge any updates since you forked into fix or feature branch, then push
   your fix or feature branch to your GitHub, and submit a pull request to
   the `develop` branch on the BEHR repository.

## Testing

As the BEHR devs review your pull request, they may ask for evidence that the modification
works as expected. The best way to test that the change does not introduce an unexpected
changes is to run `unit_test_driver` in `BEHR_main/Production tests`, which will produce
BEHR data from the starting OMI NO2 and MODIS data for certain days that test

## Style guide

Contributors should be familiar with the following conventions and follow them as 
much as possible in any new code:
