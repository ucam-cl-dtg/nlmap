AUTOMAKE_OPTIONS = gnu
lib_LTLIBRARIES = libNLMAP.la

libNLMAP_la_SOURCES = 	AngulationFunction.cc \
			MagnitudeSorter.cc \
			GaussJordanEliminator.cc \
			LaterationFunction.cc \
			LaterationSorter.cc \
			LaterationData.cc \
			NonLinearModel.cc \
			IterativeModeler.cc \
			MultiAngulation.cc \
			MultiLateration.cc \
			FitFunction.cc \
			FitData.cc



nobase_include_HEADERS = 	NLMAP/AngulationFunction.hh \
				NLMAP/MagnitudeSorter.hh \
				NLMAP/GaussJordanEliminator.hh \
				NLMAP/LaterationFunction.hh \
				NLMAP/LaterationSorter.hh \
				NLMAP/NLMAPExceptions.hh \
				NLMAP/NonLinearModel.hh \
				NLMAP/Parameters.hh \
				NLMAP/IterativeModeler.hh \
				NLMAP/MultiAngulation.hh \
				NLMAP/MultiLateration.hh \
				NLMAP/Positioner.hh

libNLMAP_la_LDFLAGS = -version-info 0:0:0

bin_PROGRAMS = example
example_SOURCES = main.cc
example_LDADD = libNLMAP.la
noinst_PROGRAMS = example
