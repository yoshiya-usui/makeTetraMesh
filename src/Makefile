CXX           = icpc
CC            = icpc
CXXFLAGS      = -O3 \
                -qopenmp \
                -D_USE_OMP \
                -D_THICKNESS_FUNC \
                -DNDEBUG 
DEST          = ./
LIBS          = -mkl=parallel
OBJS          = AnalysisDomain.o \
                BoundaryCurve.o \
                BoundaryCurveInner.o \
                BoundaryCurveList.o \
                BoundaryCurveOuter.o \
                BoundaryCurveSubInner.o \
                CoastLine.o \
                CoastLineList.o \
                Control.o \
                Ellipsoids.o \
                main.o \
                LakeList.o \
                Node.o \
                NodeList.o \
                ObservationLine.o \
                ObservationPoint.o \
                ObservingSiteList.o \
                OutputFiles.o \
                TopographyData.o \
                TopographyDataList.o \
                Triangle.o \
                TriangleList.o \
                Util.o
PROGRAM       = makeTetraMesh

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(PROGRAM)

clean:;		rm -f *.o *~ $(PROGRAM)
