
for sample in range(0, 10000):
    sample_dir = "sample_" + str(sample)
    os.system( 'rm finaliSSHadrons.dat' )
    os.system( 'mkdir ' + sample_dir)
    #run sampler
    os.system( './iSSSamplingOnly' )
    #copy results to sample dir
    os.system( 'cp finalliSSHadrons.dat ' + sample_dir + '/.')
