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
			FitData.cc \
			RANSAC.cc \
			TriAve.cc \
			MinTriLat.cc \
			MaxTriLat.cc \
			LinearisedLS.cc \
			RobustLS.cc \
			LeastSquares.cc



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
				NLMAP/RANSAC.hh \
				NLMAP/TriAve.hh \
				NLMAP/MinTriLat.hh \
				NLMAP/MaxTriLat.hh \
				NLMAP/LinearisedLS.hh \
				NLMAP/RobustLS.hh \
				NLMAP/LeastSquares.hh

libNLMAP_la_LDFLAGS = -version-info 0:0:0

bin_PROGRAMS = example ransac_test triave_test ls_test robust_test

example_SOURCES = main.cc
example_LDADD = libNLMAP.la

ransac_test_SOURCES = ransac_test.cc
ransac_test_LDADD = libNLMAP.la

triave_test_SOURCES = triave_test.cc
triave_test_LDADD = libNLMAP.la

ls_test_SOURCES = ls_test.cc
ls_test_LDADD = libNLMAP.la

robust_test_SOURCES = robust_test.cc
robust_test_LDADD = libNLMAP.la

noinst_PROGRAMS = example ransac_test triave_test ls_test
