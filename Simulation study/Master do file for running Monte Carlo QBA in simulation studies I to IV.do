/***********************************************************************************************************************************
                              MASTER DO-FILE FOR RUNNING SIMULATION STUDIES I TO IV FOR MONTE CARLO QBA

RUNS DO FILES FROM FOLDER "Do Files" AND POST RESULTS TO A FOLDER CALLED "Results". 
							  
***********************************************************************************************************************************/
version 18

* SET WORKING DIRECTORY; E.G., 
* cd "C:\Simulation study\"

/************************
   SIMULATION STUDY I
*************************/
/************************************************************
  SCENARIO A: Y BINARY, X BINARY, U=(U1) CONTINUOUS
	VARY LEVEL OF INFORMITVENESS OF BIAS PARMETERS' PRIORS
	VARY PRIOR USED FOR RESIDUAL SD eta
	VARY NUMBER OF MONTE CARLO REPLICATIONS 
*************************************************************/
* VERY INFORMATIVE; UNIFORM DISTRIBUTION FOR BIAS PARAMETER eta; REPEAT FOR 100, 400, AND 10,000 MONTE CARLO REPLICATIONS 
set seed 1058
local samplesize 1000
local numsimdatasets 500
local priorlevel 3
local uniform 1
foreach numMCsteps of numlist 100 400 10000 {
	run "Simulation study I\Scenario A\Monte Carlo\Do files\Scenario A - logit - runs sim study I for Monte Carlo QBA.do" /// 
	`samplesize' `numsimdatasets' `numMCsteps' `priorlevel' `uniform'
}

* INFORMATIVE; UNIFORM DISTRIBUTION FOR BIAS PARAMETER eta; REPEAT FOR 100 AND 400 MONTE CARLO REPLICATIONS 
set seed 1115
local samplesize 1000
local numsimdatasets 500
local priorlevel 2
local uniform 1
foreach numMCsteps of numlist 100 400 {
	run "Simulation study I\Scenario A\Monte Carlo\Do files\Scenario A - logit - runs sim study I for Monte Carlo QBA.do" /// 
	`samplesize' `numsimdatasets' `numMCsteps' `priorlevel' `uniform'
}

* INFORMATIVE; UNIFORM DISTRIBUTION FOR BIAS PARAMETER eta; 10,000 MONTE CARLO REPLICATIONS 
set seed 1841
local samplesize 1000
local numsimdatasets 500
local priorlevel 2
local uniform 1
local numMCsteps 10000
	run "Simulation study I\Scenario A\Monte Carlo\Do files\Scenario A - logit - runs sim study I for Monte Carlo QBA.do" /// 
	`samplesize' `numsimdatasets' `numMCsteps' `priorlevel' `uniform'
																		   
* INFORMATIVE; GAMMA DISTRIBUTION FOR BIAS PARAMETER PRECISION 1/etasq; REPEAT FOR 100, 400, AND 10,000 MONTE CARLO REPLICATIONS
set seed 2102
local samplesize 1000
local numsimdatasets 500
local priorlevel 2
local uniform 0
foreach numMCsteps of numlist 100 400 10000 {
	run "Simulation study I\Scenario A\Monte Carlo\Do files\Scenario A - logit - runs sim study I for Monte Carlo QBA.do" /// 
	`samplesize' `numsimdatasets' `numMCsteps' `priorlevel' `uniform'
}

/**************************************************************
  SCENARIO B: Y CONTINUOUS, X CONTINUOUS, U=(U1,U2) CONTINUOUS
	VARY LEVEL OF INFORMITVENESS OF BIAS PARMETERS' PRIORS
	VARY PRIOR USED FOR RESIDUAL SD eta
	VARY NUMBER OF MONTE CARLO REPLICATIONS 
*************************************************************/
* VERY INFORMATIVE; UNIFORM DISTRIBUTION FOR BIAS PARAMETER eta; REPEAT FOR 100, 400, AND 10,000 MONTE CARLO REPLICATIONS 
set seed 1431
local samplesize 1000
local numsimdatasets 500
local priorlevel 3
local uniform 1
foreach numMCreps of numlist 100 400 10000 {
	run "Simulation study I\Scenario B\Monte Carlo\Do files\Scenario B - regress - runs sim study I for Monte Carlo QBA.do" ///
															`samplesize' `numsimdatasets' ///
															`numMCreps' `priorlevel' `uniform' 
}

* INFORMATIVE; UNIFORM DISTRIBUTION FOR BIAS PARAMETER eta; REPEAT FOR 100, 400, AND 10,000 MONTE CARLO REPLICATIONS 
set seed 250319
local samplesize 1000
local numsimdatasets 500
local priorlevel 2
local uniform 1
foreach numMCreps of numlist 100 400 10000 {
	run "Simulation study I\Scenario B\Monte Carlo\Do files\Scenario B - regress - runs sim study I for Monte Carlo QBA.do" ///
															`samplesize' `numsimdatasets' ///
															`numMCreps' `priorlevel' `uniform' 
}

* INFORMATIVE; GAMMA DISTRIBUTION FOR BIAS PARAMETER PRECISION 1/etasq; REPEAT FOR 100, 400, AND 10,000 MONTE CARLO REPLICATIONS
set seed 4341
local samplesize 1000
local numsimdatasets 500
local priorlevel 2
local uniform 0
foreach numMCreps of numlist 100 400 10000 {
	run "Simulation study I\Scenario B\Monte Carlo\Do files\Scenario B - regress - runs sim study I for Monte Carlo QBA.do" ///
															`samplesize' `numsimdatasets' ///
															`numMCreps' `priorlevel' `uniform' 
}

/************************
   SIMULATION STUDY II
*************************/
/**************************************************************************
  SCENARIO A: Y BINARY, X BINARY, U=(U1) CONTINUOUS
	VARY LEVEL OF STRENGTH OF ASSOCIATIONS BETWEEN C AND U 
	VERY INFORMATIVE PRIORS; UNIFORM DISTRIBUTION FOR BIAS PARAMETER eta
	100 MONTE CARLO REPLICATIONS 
***************************************************************************/
set seed 1558
local samplesize 1000
local numsimdatasets 500
local priorlevel 3
local uniform 1
local numMCreps 100
foreach strength_UC of numlist 0 2 {
	run "Simulation study II\Do files\Scenario A - logit - runs sim study II for Monte Carlo QBA.do" ///
																								   `samplesize' `numsimdatasets' `numMCreps' ///
																								   `priorlevel' `uniform' `strength_UC'
}

/**************************************************************************
  SCENARIO B: Y CONTINUOUS, X CONTINUOUS, U=(U1,U2) CONTINUOUS
	VARY LEVEL OF STRENGTH OF ASSOCIATIONS BETWEEN C AND U 
	VERY INFORMATIVE PRIORS; UNIFORM DISTRIBUTION FOR BIAS PARAMETER eta
	100 MONTE CARLO REPLICATIONS 
***************************************************************************/
set seed 1735
local samplesize 1000
local numsimdatasets 500
local priorlevel 3
local uniform 1
local numMCreps 100
foreach strength_UC of numlist 0 2 {
	run "Simulation study II\Do files\Scenario B - regress - runs sim study II for Monte Carlo QBA.do" ///
																								   `samplesize' `numsimdatasets' `numMCreps' ///
  																						           `priorlevel' `uniform' `strength_UC'
}

/************************
   SIMULATION STUDY III
*************************/
/**************************************************************************
  SCENARIO C: Y BINARY, X BINARY, U=(U1) BINARY
	VARY LEVEL OF STRENGTH OF ASSOCIATIONS BETWEEN C AND U 
	PREVALENCE OF U AT ~20%
	VERY INFORMATIVE PRIORS; UNIFORM DISTRIBUTION FOR BIAS PARAMETER pi
	100 MONTE CARLO REPLICATIONS 
***************************************************************************/
* ASSOCIATIONS BETWEEN U AND C SET AT VALUES OBSERVED IN REAL DATASET
set seed 1412
local samplesize 1000
local numsimdatasets 500
local priorlevel 3
local uniform 1
local numMCreps 100
local strength_UC 1
local pilevel 1
run "Simulation study III\Do files\Scenario C - logit - runs sim study III for Monte Carlo QBA.do" ///
                                                                                 `samplesize' `numsimdatasets' `numMCreps' ///
                    															 `priorlevel' `uniform' `strength_UC' `pilevel'

* INDEPENDENCE BETWEEN U AND C 
set rngstate XAA54602691f086380ce3b62b2b9ea3d3b3598589abc4b2a0783a0a3dd981be88777e4a140038bb9f2dfc8a725e78124f284908917eb8d6dd0ff35d0bc8cea1e6af4186219c6ce481ac7c5dad8e661a1967ea30db99f59262722e5d704a93aadfb698c7864bf91b56a5f7aa55ec38d249d00f8c4888ab7518c23bba1c41598a08660ac8b68061e1a09eb1e3e24eaefb58bdb3cc2767199cbff13ee8774b9c0e03b7af05377247428225feca8993bea954a7f2409710074fe1de96def83f9e3523675bb7e286e7b4816b55e0c9ff2e5dcfa7b8d783aa6ca5676d6947db3af2ad4a2ae7f8b359720ec7f45729df8f35f76c3d438a35f786bcc2e3701edab0eeb0d7aacbdf072a18ecf2cbc805f3c72aae20a6b3a10257f60da3bd4211e54263ec500afe33abb4aa4d4c4b219ef2bdbb3d27c133affcd31aafddbb1e06964f3d0464220433437696888d1ae58ab223db73c4fc402e3ba8060a0aede57b3008ffba5e588a860ab6ea46019f43c94f92b854c24949a779a2e978b801bbff9b1cce6d174b54279128081d0623d57aacaa3f6ad75bef5c822f201b677a76fa4367366bc67977eef70c1557f14408c53dbd429d3f2665f12dce1d8d9df86626d994cab6cfe65e1a75aa69582117680e623a71f1f5b2e44ba01b6eb2260c61d8ea4950ff734ea9080096bbda39fac84dc6cfe869400fd825d0eb75427a2bea05dc32f7af3b649cc227ee05e05dcf42ff54bd9a936e7ab301900b9240b5acc84f6c2e17372c399eecf44adcd24c2d997674e244775244bd57c4aabad7e73e5df3989e91515acc9591119fa899dca98a12a2e78bea87a7c3998ce59d6c59fe5dd7ea0b57bf6678e1d184a17ff2f5b3e890d2cdd67ba2b7a29e3e61f7111eb27f352f92d2f2c2b9a4a743e336b3a6ffc03cac4afa559ad05d4aa233c13e5bec6fa97ef89e6193703640514fa96008197d948da99fd49fe19b63d66a4bee5e52e4bf154bf4b2390d55cca2a9c0d9d660424777a9ad674cb6b0933cf3ea7bf0f31f70e671b4f6b329dbb91478f0372fd8fb3142d9260461763ca66f1ef7274637e73312e9244d827760df590f4bca77e06affd6cfbf92d55ea08446bec0ed6e752419559b348bf477733ea0d840b611d025b82bef49e343b144f3afa879dd7f2e28dc4041dab0fb3daa94ab8631fdf019c6d6b5479048cafcd156f6aa7897b1dc0210b10591669e788569c4f9ac73171d85d542015c482a3d014caa8fbcb480ac8d635e62ef0e90be10862395983f8e3f0d9b2a86c1d25c0cb764c3dfcca0f901ed4eb73b2ef14a5ece1a04ee583d07ec50f1856c9b4fd9647d8a3725843b6d73329350b576cfb407819dbe34855b12c76b69f93587a6d36e78a17400900f01dbd998cb323e0e62b9136855c92ab0c6ae7821a7d809f1fab70c16e9823b60493ffcc95c07d6d0718c933b42b107ea24c4b3d43281ce2bdcbdce2e07b100ce97ad2ee7b3406cf2ab169bf9f3a9b8df2b23dffc122170c7ba5863496cb2f835382662e2ef809a5f4c264fc308b5583227c4159b793b8e2d90de89022915c5eb97fb7ea6e7d9b11fb133209ce69aba93a1da587a1b5b139fd42d0caf33d4ac57a22a97560d255d886cd41e55094a3f4c58cc0bc234c36717b85d918f1d84a51ca8fc399176fc793e8f1fb47d30b79dbed213b4ca66450154e138472f3403e87046100bf0515c9407e3e2900df8a1ccbb6cbee88b66491fd6806df0c1702b1fd6c5814b8be50cd55fe187963a615b1eb372122d454be82c94edd987efc1c69b2ad8f5da7e53b1a010ce84bccb4eaf01e4fde1e0553744310ed1f2d19d987a0ffd414a71f3a500811774aa0ee507e913c56fc99c31ae8bc6ab69fa8a7f666061aa53e57fa5366b541ab6d75c9fc51cef6692205d89ce78aab7ade9cbb6c38284eacc6eb380a71a3310268e1c0f41a7ade9f8b7e59366eb3ceaa078d2620aa65739a84a8bb5bb3fd47be336813e5cc1d412a297ac182d8c060168813c66f6ca23e3ce9c98976666e25aab7cdaf65f333139677bec2b5abda0bd409d76d29bfa6f397aea6efa7448164375d58fb09eac55b6dfdfa2641c16869c69e830dc5344bcd69e9cad60d7d50de02aef2583315ab4a72c305782fc16ec1a622719f54e549d5fac4e72c14fcef0565420e11b8a2074dbae6247da8dae15bd153ecd281c795ecb0273cf60070e31c603c1b86e485317af8a24a9f0312d9d54f4dc685a18cd912c31ba88fa59c4964f21a4418b490c677b690e5ead32211bfbba0b4bbbb98d3103cbd96877793a2e1d7dd77feace7a975a7d68457453a933a19811c45be1935da0385e835d818162cabdd45af4763ffc2b422a1ab98cff4daeff1f0d931d0f4ebaddfb7dc54c76c4aa352a2b0f4441b585d337a4f5ccc25c6c92174ea429877203b04e61f748e3e66c22136227bc834b6251e18b6718469894177b74c0cfedb4c56557522f0507da3e28bc612e089226ac0be55d665aa7eb44897124355cf7226912a073d1fa6cc205a9ecd440e4733a070defe39577d760440d48c2fcc1125ef47cc348393defde762923a313cbda9a6876478d8de4943899841aefaffd0dde6d686d630ff866c0a16e03f41ff1f815052a206392d60311e084c6ca2e588c4504176e022e410d0bad9c3107f78ab4b14668ee5edc380d759d32e0366b20bf275f3a4cace251a907d073a12ea36e55999ef89773c2cd455e790d0babc50ad0eb5492c8d5ecb647c2c5ce1e5fd42cd5f783d71bf160902cb5c39126a719c44d96d4523056c06492c4a08398cec0a22e569e2f27d488e5e21b523bab2bc5bbd029f0ed7df928cf5a8c1b580b4bb8b2a2028209c568c20dc80262d061d57442ea4390131094817d3fd1c994646e1df79ac2b76440debf70c559992035c28b046575fcf1243137afde42698c529ff2531278f31f80c0166b3ff2c7a2a5f0e9a45389d3fcd37a774ea00c10d7578b01bdd559f3a8c18ee9dde8076e534a5925643c75f20df4c056f833e05a37c4407f37a926c170e20bfcb19442b5fdb8e5bdc20479d67756e18c2530fbb3098cdad7af154e5e54447b867e7452a956dfeb64cbd5f35a6523c044a3c59aa2a194a1679f61a554790c481c30590ebf180d167b23b0ea82897c978c8059a32bf3cc790dbf002dafd57175a134b4056345a2005544703d6282df616b22559e17852835c42b8848dc69c308ae9f11046261b9162ff767b273d432863fb5fbae55ed9742bfe5a5ef479be83e3719399bf28bc6f40eb2a1060d76435a7798e8da61841b75b869da04529005150e556e49678c94e2300de7d618d0acd6f08c0d82af87755f5757a8deacb1fc470fba110f6879a9518149d8b1eebffcfef2b3e0c3652b71be055de952e3f2ae6b49f21b8193251e6a4da1661f3dd9dd9e7f4e6b6f1ab953aa5690468614ce6ea6bbbf32826ceb05a67ea0c4613e247b6adee9f62fdeef926cbacaa820853bea7f3edb992123cf7e4a497ac1f2410dd84110c5086d1cfd810001000000ee4c33
																				 
local samplesize 1000
local numsimdatasets 500
local priorlevel 3
local uniform 1
local numMCreps 100
local strength_UC 0
local pilevel 1
run "Simulation study III\Do files\Scenario C - logit - runs sim study III for Monte Carlo QBA.do" ///
                                                                                 `samplesize' `numsimdatasets' `numMCreps' ///
                    															 `priorlevel' `uniform' `strength_UC' `pilevel'

* ASSOCIATIONS BETWEEN U AND C DOUBLE THE VALUES OBSERVED IN REAL DATASET
set rngstate XAA832d6c1dd71bbacaad7ca3d83518a9724107fd9aafd1d92ae57116c4599714beeb4139e23f0e99fc724deddf49cde745168825747acf524965f35bdb6c9574edf6ae0d8b7c1a39b13cd6a0f12d56d56333ce1e4889cd572e35987c3f43e6327a699e33c3643e4acaae91edc953a88ce1016b2ce596b30a23c210480944c9f68ef8e2d01a590b37b4707893c18ee15c42e62aa4bb7a8552c07ffb9c961a67e4b44178286ad3a6882af56e44738c3b4cd51b5c4062409ef9ae39a9d25c08176665a7fb19e0ad9a167b9335f01a47c65d70517493f978de1296c5037431eb463357cfa82a6b21ca71db90caa7b114372448c8ea551e52b8fbf05dc23385305c6155402c14ae3eaed903dd55a0127461a945c0e9d4305df4634aff5bc4a9862bfa8f8e88550b4fcbfcba7cc10ec15dce5b236cc689301be8e02194782b6f1f4200af71fb72ea62cc5973194391949b310bb96760864bb3e51bb6ffb6c1ae9787734dbbb617cbebf407c7f1c901b37970059f59a98e5f4d6ae6b54910192b8927eb162bb2a75dd7101f249b13c0fec4a96fcb115b6e06fc842ebaab72a9f221cd0a0392d23c62114e2e32768d7756adabf94a3274eee6d9262dd25231cff7ab8ab50926b4e427929597ee6aecd840e96e467be8af57f8b9717d15194c7424dff3928b99bec642fbe7fdd981ac594ccbb57723a916d52f55874f568f69b99fa541ca35abc12f831a6181e726bc866798bc03a2b17b705101a9e7d68f03b43d00129d3cd025e406b071b38875e1b7b2fe1e7cbacdb780b4f434ac7f69716fa9af6b094c0510ce6ab4a1ca5e916df344f03086b233b2807e6d0f45933bba79847b9b3a29001e0bb5ea7b3bfe25029a8620d268f72fb8e615e9bd304617e39df423aac5c7855414866e790f9f7f1cf9ea1369c516f157b705baef06f4e5c09526e25a76268715337fc138ccadf8e80cb3e2ea9402e8f34d7933b963ec3b838f079f6a46ccb4076f5f6a5bd2e36c8119f8d0ab27ad35e7a22a441f32e53ed2861ba11cb83fc8f190ba541b78b4137ad8c6d6ebf5915d3e32215d58af232c9b6b85bda2aa6a4470170011900e45e70f7d8e84d870bc52fe6b21426b4af9b45078ba8f1082383554c6d059baa25fc9a36c5c17ef4f01e7c60d145ffc6778a8de0bbee6c1bf6e1011a07e8d82826d7668fc97b3f8043beab8e921322df6c785e3efb53ed57b89b793223b7e89d6c720c5e87e089ccb4d9e439a07cda8921a738e86c7d37db2af8e10b7181f37b24594987ddcd7171a7aeaaba7808bf2b6742cb8846e8b323f133ee3b6137e3a9c855b0b7f560185f250ec39a00bee324b0f7cc38aec9d8d81c82a373d3091a3d637ee7e4e09f1c4a37c935f17c8bab35a632c700f5835f97b96ba6356227639ac55f91762a1e98ed51c282594d65bc76511e6da7df8bbda0847251f610cb896b6c1a5e6708018bb5131b26f435e09089a08283b4401c5e8fb67189293ab9850a60d246a2237c767e818950feadc7fcd6cc3f34926743ba2304fd329a76748d5c11153d45ccd2367a6fa35fa75c793a6eade0901334ed4e3d198428dabca6d7f58c3dd8033ea9f0088e5267d8532653e200299aae0faf9a09d58a332cf0f913ceb69322544e5ece29b581d529e2df6c89876f1ba05afb127be6d1a331a1a063bd56cc2048bb6f306ac80af0c1d4c0ea496618c8e2ba92342cc065fd6f17b272162db64f90098a0ad77fafe6e9acdfcca45d7416eeef4031f95f3007a6ff8ffdc6992326d1b0feae617b13ff70bd80653f6ceadc81fcfeac1129e614c8d7fa50aa0d1a27156c1630bde515680bcd03d087a3c729294d18352f2998e48fdb725928e2660a99a2a5d297db738e3047b5217fb01e234dad15f42b0d914609fd24f8a7af1ca0ce2ea89fbb2ab42da5d139fe9709e57de190e04e321e124051380dba4d514b0056d7af7274b64a10b7018b0501cd698e5e6788a8cd392187101485dd2a579ea51245aa7b61e522b083764692a68f9fa0b8ddcb01928137d7c3ce7207834c931e95ca794849f7d343f737b9815c83214476bac46b07b5f0397194c19d19faecb778d1ad7029ac2434c39d93d5bf0e17e1f8bb608cb5959c19650ef04c5f36a0569c40b53f064667e42247c6ad9a44d79a50c2c9038025e55a4a2eda7e4f4d54db9daf4f05d61d738886652479558277fa44459dd64ca4f70b2a2bed2d51d1fd2520b878d0367cb34d5b638148afc18514fa2b5ba5114b5d6ad4e1aedc6cddfc1234d7982dc7a4a4a7c20b3c15ea105fa495add0bb57f5518aa4867db433af2565d816213d085bd2e00dbd596789953f06e231a1827cb7428ebd37a28c7c1f40cd99b770f1e29d9aee4090d594b832f7547f4613e87d9d07e093b5443c78ce050624a7e2912dd7272f1f287b3b911be03bebe648d63b900be95dfd4d89df20efaaa1d26431bd33543f2318874e5a4af162cea0bae0db5a9677db089030e36bb36285931331d965c810fc8a9f56681557f365193f4380b4a29a45356a54f0bf7f8e9f6a4c2245ed9d05d58381bf590659c5e45a195d8d5220d275ab2760fb544edad0d313656e625fe38b03bbbe85fd0ded4bec65926344a27f57029de639293a6afa704ebc8b0181a4af5a7a971c78f69bf5524a85ce7a12587f782fc9c77539d7e7f04605bae95a5db6a58bf9efe4da8331eee71444042c79ccd6835ab9d31296f92a91b5520db088a2fd53fdb97ad4cf6d0f91e2cbfa68773d86f035a6cd005908bfb50d628f6595ecd824d72deadc960d21343d986a5e0428c9d9a5009115cc5204b004b6fa26e1387a240ee0aff3d99c5b86029065b6cd1e7e233c2158ea67797e3044503a9e0695544b0ccea9cdd69947bfba4f30719cc246114e899e63a698723011f2af8811b68711a0972a406fecf11bded0b70e1ccd2d14d0afdb12d400cfcc8fe793cb3aa7fc380abcad08dfe9549064f8987ed6bfe93d9cce4b0aeb5aefd52384c94e04c1764ca796fb9f14e5a30d47a2d54f8347fd594f935a164b45da1ab2425dac57bb7e05035d6f842e0521a9776639f4475e48ba3686ad08a9711a8dc14464cc69b5c0655e8b2af85c481ada7ebfcfc811601c9bcf90e6ea75628c74387e4d8886977ec619c1bd6a5ac8ab2ce89c8dcfba72cc12b58152e86110444ef2497edd16505f2ad36f2142f6709baf5ebdd626f5c84e6105fb6964ac98dde91a17252de1ea26feb0c64f86f815d362495a575619d0295ac9cc18ef74b78901ff1a3eaffbb676a0a80eda963731b962ae8b345db29a62a2b3f63feb2206d44986b57e79036238848775d870a1963b31126394b4a52c02eb1c8bf3ebfe944292f32847f3499fc9e81c599b7acd37d0825989589b1b471232d4ea75316cc57b90c21be6c6ba2587263a247d4a41eec3b6493f86a78374d13af6a2ac6658e15ed6f74496bbdd822b60ecb62b3fc6809058d89d8611cb55298eabf03bc7674a8b15e3b8f82a031515cd4652e4b6ddf4f90a892489600010000003e0148
local samplesize 1000
local numsimdatasets 500
local priorlevel 3
local uniform 1
local numMCreps 100
local strength_UC 2
local pilevel 1
run "Simulation study III\Do files\Scenario C - logit - runs sim study III for Monte Carlo QBA.do" ///
                                                                                 `samplesize' `numsimdatasets' `numMCreps' ///
                    															 `priorlevel' `uniform' `strength_UC' `pilevel'
																				 
/***********************************************************************************
  SCENARIO D: Y CONTINUOUS, X CONTINUOUS, U=(U1,U2) BINARY
	VARY LEVEL OF STRENGTH OF ASSOCIATIONS BETWEEN C AND U, AND BETWE 
	Pr(U1=1)~14.1% AND Pr(U2=1)=18%
	VERY INFORMATIVE PRIORS; UNIFORM DISTRIBUTION FOR BIAS PARAMETERS pi1 AND pi2
	100 MONTE CARLO REPLICATIONS 
***********************************************************************************/
set seed 1955
local samplesize 1000
local numsimdatasets 500
local priorlevel 3
local uniform 1
local numMCreps 100
local pilevel 1

foreach strength_UC of numlist 0 1 2 {
run "Simulation study III\Do files\Scenario D - regress - runs sim study III for Monte Carlo QBA.do" ///
																		  `samplesize' `numsimdatasets' `numMCreps' ///
                                                                          `priorlevel' `uniform' `strength_UC' `pilevel'
}

/************************
   SIMULATION STUDY IV
*************************/
/**************************************************************************
  SCENARIO E: Y CONTINUOUS, X NOMINAL, U=(CONTINUOUS U1, BINARY U2)
	STRENGTH OF ASSOCIATIONS BETWEEN C AND U AND U1 AND U2 SET TO OBSERVED VALUES
	VERY INFORMATIVE PRIORS; UNIFORM DISTRIBUTION FOR BIAS PARAMETERS eta1 AND pi2
	100 MONTE CARLO REPLICATIONS 
***************************************************************************/
set seed 1015
local samplesize 1000
local numsimdatasets 500
local priorlevel 3
local uniform 1
local numMCreps 100
local strength_UC 1
local pilevel 1

run "Simulation study IV\Do files\Scenario E - regress - runs sim study IV for Monte Carlo QBA.do" /// 
																		  `samplesize' `numsimdatasets' `numMCreps' ///
                                                                          `priorlevel' `uniform' `strength_UC' `pilevel'

/**************************************************************************
  SCENARIO F: Y NOMINAL, X CONTINUOUS, U=(CONTINUOUS U1)
	STRENGTH OF ASSOCIATIONS BETWEEN C AND U 
	VERY INFORMATIVE PRIORS; UNIFORM DISTRIBUTION FOR BIAS PARAMETERS eta1
	100 MONTE CARLO REPLICATIONS 
***************************************************************************/
set seed 250613
local samplesize 1000
local numsimdatasets 500
local priorlevel 3
local uniform 1
local numMCreps 100
local strength_UC 1
run "Simulation study IV\Do files\Scenario F - mlogit - runs sim study IV for Monte Carlo QBA.do" /// 
																		  `samplesize' `numsimdatasets' `numMCreps' ///
                                                                          `priorlevel' `uniform' `strength_UC' 

/**************************************************************************
  SCENARIO G: Y SURVIVAL, X CONTINUOUS, U=(CONTINUOUS U1)
	STRENGTH OF ASSOCIATIONS BETWEEN C AND U
	VARY LEVEL OF INFORMITVENESS OF BIAS PARMETERS' PRIORS
	100 MONTE CARLO REPLICATIONS 
***************************************************************************/
* VERY INFORMATIVE
set seed 1732
local samplesize 1000
local numsimdatasets 500
local priorlevel 3
local uniform 1
local numMCreps 100
run "Simulation study IV\Do files\Scenario G - Cox PH - runs sim study IV for Monte Carlo QBA.do" /// 
																		  `samplesize' `numsimdatasets' `numMCreps' ///
                                                                          `priorlevel' `uniform'  
* LESS INFORMATIVE
set seed 903
local samplesize 1000
local numsimdatasets 500
local priorlevel 2
local uniform 1
local numMCreps 100
run "Simulation study IV\Do files\Scenario G - Cox PH - runs sim study IV for Monte Carlo QBA.do" /// 
																		  `samplesize' `numsimdatasets' `numMCreps' ///
                                                                          `priorlevel' `uniform'  
																 
* END OF DO-FILE