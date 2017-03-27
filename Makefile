all dep clean indent tests::
	cd gtest/lib && make $@ && cd .. \\
	cd testlocalopt && make $@ && cd .. \\
	cd testsearch && make $@ && cd .. \\
	

