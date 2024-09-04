



default: TestRun

TestRun: Execute
	wget https://cernbox.cern.ch/remote.php/dav/public-files/7ST6Uclh5VjJnt7/workshopexercise10k_1.hepmc
	./Execute --Input workshopexercise10k_1.hepmc --Output test.out
	# cat test.out | TextToTree test.root 6 "JetPT:JetEta:JetPhi:CAELoss:TimeELoss:KTELoss"

Execute: EvaluateEnergyLoss.cpp CATree.o TauHelperFunctions3.o
	g++ EvaluateEnergyLoss.cpp -o Execute TauHelperFunctions3.o CATree.o -std=c++17 \
		`fastjet-config --libs --cxxflags`

CATree.o: CATree.cpp CATree.h
	g++ CATree.cpp -c -o CATree.o -std=c++17

TauHelperFunctions3.o: TauHelperFunctions3.cpp TauHelperFunctions3.h
	g++ TauHelperFunctions3.cpp -c -o TauHelperFunctions3.o -std=c++17
