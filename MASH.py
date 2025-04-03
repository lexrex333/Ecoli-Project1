def MASH(one,two):
	mash_command = ["mash dist ", one, two]
	print(mash_command)
	#subprocess.run(mash_command)
test = "Test"
test_two = "Test_two"
MASH(test, test_two)
