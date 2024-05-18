with (import <nixpkgs> {});
let
	my-python-packages = python-packages: with python-packages; [
		pandas
		numpy
		scipy
		matplotlib
		toml
		geographiclib
		cartopy
		eccodes
		obspy
		setuptools
	];
	python-with-my-packages = python310.withPackages my-python-packages;
in
mkShell {
	buildInputs = [
		python-with-my-packages
		gfortran8
		rocmPackages.llvm.clang
		fftw
		gnumake
	];
}
