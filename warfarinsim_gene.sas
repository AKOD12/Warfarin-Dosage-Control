options nonotes;
options nodate formdlim="-";

data temp1;
	/****Static Variables****/
%Let nsim=10000;
	*number of simulations;
	%Let nstudy=50000;
	*study population;
	%let percentAMIO=0.0588;
	*percent of people taking Height;
	
	*African-Americans;
	%Let meanAAage=6.2;
	*mean age of african-americans in decades;
	%Let sdAAage=0.9;
	*standard deviation of african-american age in decades;
	%Let percentAA=0.425;
	*percentage of african-americans in study;
	%Let meanAAweight=90.6;
	*mean african-american weight;
	%Let sdAAweight=22.9;
	*standard deviation of african-american weight;
	%Let meanAAheight=170.9;
	*mean height of african-americans in centimeters;
	%Let sdAAheight=10.5;
	*standard deviation of african-americans in height;
	
	*Europeans;
	%Let meanEUROage=7.8;
	*mean age of europeans;
	%Let sdEUROage=0.10;
	*standard deviation of european age in decades;
	%Let percentEURO=0.575;
	*percentage of europeans in study;
	%Let meanEUROweight=79;
	*mean european weight;
	%Let sdEUROweight=20.5;
	*standard deviation of european weight;
	%Let meanEUROheight=173;
	*mean height of europeans in centimeters;
	%Let sdEUROheight=7.5;
	*standard deviation of europeans in height;
	
	*Afrcian-American CYP2C9 Diplotype Frequencies;
	%Let AAprob6homo=0.0001;
	%Let AAprob3homo=0.0001;
	%Let AAprob5homo=0.0002;
	%Let AAprob3and6=0.0002;
	%Let AAprob11homo=0.0002;
	%Let AAprob5and6=0.0002;
	%Let AAprob11and6=0.0002;
	%Let AAprob3and5=0.0003;
	%Let AAprob11and3=0.0003;
	%Let AAprob2and6=0.0004;
	%Let AAprob11and5=0.0004;
	%Let AAprob2homo=0.0005;
	%Let AAprob2and3=0.0005;
	%Let AAprob2and5=0.0006;
	%Let AAprob11and2=0.0006;
	%Let AAprob6and8=0.0010;
	%Let AAprob3and8=0.0016;
	%Let AAprob5and8=0.0017;
	%Let AAprob11and8=0.0018;
	%Let AAprob2and8=0.0031;
	%Let AAprob8homo=0.0044;
	%Let AAprob1and6=0.0133;
	%Let AAprob1and3=0.0202;
	%Let AAprob1and5=0.0222;
	%Let AAprob1and11=0.0240;
	%Let AAprob1and2=0.0398;
	%Let AAprob1and8=0.1151;
	%Let AAprob1homo=0.7469;
	
	*European CYP2C9 Diplotype Frequencies;
	%Let EUROprob6homo=0.0000;
	%Let EUROprob5homo=0.0000;
	%Let EUROprob3and6=0.0000;
	%Let EUROprob11homo=0.0000;
	%Let EUROprob5and6=0.0000;
	%Let EUROprob11and6=0.0000;
	%Let EUROprob3and5=0.0000;
	%Let EUROprob2and6=0.0000;
	%Let EUROprob11and5=0.0000;
	%Let EUROprob2and5=0.0000;
	%Let EUROprob6and8=0.0000;
	%Let EUROprob5and8=0.0000;
	%Let EUROprob11and8=0.0000;
	%Let EUROprob8homo=0.0000;
	%Let EUROprob1and6=0.0000;
	%Let EUROprob1and5=0.0000;
	%Let EUROprob3and8=0.0002;
	%Let EUROprob11and3=0.0002;
	%Let EUROprob2and8=0.0004;
	%Let EUROprob11and2=0.0004;
	%Let EUROprob1and8=0.0023;
	%Let EUROprob1and11=0.0027;
	%Let EUROprob3homo=0.0050;
	%Let EUROprob2homo=0.0159;
	%Let EUROprob2and3=0.0179;
	%Let EUROprob1and3=0.1133;
	%Let EUROprob1and2=0.2017;
	%Let EUROprob1and1=0.6401;
	*African-American VKORC1 allele frequency;
	%Let AAVKORC1=0.10274;
	*Frequency of VKORC1 allele in African-Americans;
	*European VKORC1 allele frequency;
	%Let EUROVKORC1=0.41242;
	*Frequency of VKORC1 allele in Europeans;
	*African-American CYP4F2 allele frequency;
	%Let AACYP4F2=0.0772;
	*Frequency of CYP4F2 allele in African-Americans;
	*European CYP4F2 allele frequency;
	%Let EUROCYP4F2=0.29963;
	*Frequency of CYP4F2 allele in Europeans;

	%Let probhighoverdosedpgx = 0.693; *Probability of overdose if assigned high dose;
	%Let probhightruepgx = 0.307; *Probability of true dose if assigned high dose;

	%Let problowtruepgx = 0.449; *Probability of true dose if assigned low dose;
	%Let problowunderdosedpgx = 0.551; *Probability of underdose if assigned low dose;

	%Let probmedoverdosedpgx = 0.344; *Probability of overdose if assigned medium dose;
	%Let probmedtruepgx = 0.865; *0.521 + 0.344; *Probability of true (or overdose) if assigned medium dose;
	%Let probmedunderdosedpgx = 0.135; *Probability of underdose if assigned medium dose;


	/****Data Generation****/
	%let i=1;

	do j=1 to &nstudy.;
		*coding race(Race);

		if ranuni(&i+137587) le &percentAA. then
			Race=1;
		else
			Race=0;
		*1=african-american and 0=european;
		*coding age(Age);

		if Race=1 then
			Age=&meanAAage.+&sdAAage.*rannor(&i+372837283);

		if Race=0 then
			Age=&meanEUROage.+&sdEUROage.*rannor(&i+3782173917);
		*based on race, age calculations are performed;
		*coding weight(Weight);

		if Race=1 then
			Weight=&meanAAweight.+&sdAAweight.*rannor(&i+378222);

		if Weight<46 then
			Weight=46;

		if Race=0 then
			Weight=&meanEUROweight.+&sdEUROweight.*rannor(&i+758487);

		if Weight<46 then
			Weight=46;
		*coding height(Height);

		if Race=1 then
			Height=&meanAAheight.+&sdAAheight.*rannor(&i+4256651);

		if Race=0 then
			Height=&meanEUROheight.+&sdEUROheight.*rannor(&i+42987);
		*Coding African American Genotype;
		z=ranuni(&i+73927382);
		z1=ranuni(&i+328932);
		z2=ranuni(&i+10283378);
		z3=ranuni(&i+7896838);
		z4=ranuni(&i+986823);
		z5=ranuni(&i+3920937);
		z6=ranuni(&i+46782927);
		z7=ranuni(&i+789028783);
		z8=ranuni(&i+133783983);
		z9=ranuni(&i+7683297);
		z10=ranuni(&i+10283);
		z11=ranuni(&i+65387);
		z12=ranuni(&i+4532234);
		z13=ranuni(&i+654324);
		z14=ranuni(&i+76543);
		z15=ranuni(&i+56347483);
		z16=ranuni(&i+158726419);
		z17=ranuni(&i+567763290865);
		z18=ranuni(&i+9088888888843);
		z19=ranuni(&i+3752866193);
		z20=ranuni(&i+1962739020867);
		z21=ranuni(&i+729838723);
		z22=ranuni(&i+732673);
		z23=ranuni(&i+6532871679812);
		z24=ranuni(&i+5732325);
		z25=ranuni(&i+183783);
		z26=ranuni(&i+463287712);
		x1=0; 
		x2=0;
		x3=0;
		x4=0;
		x5=0;
		x6=0;
		x7=0;
		x8=0;
		x9=0;
		x10=0;
		x11=0;
		x12=0;
		x13=0;
		x14=0;
		x15=0;
		x16=0;
		x17=0;
		x18=0;
		x19=0;
		x20=0;
		x21=0;
		x22=0;
		x23=0;
		x24=0;
		x25=0;
		x26=0;
		x27=0;
		x28=0;
		
		*CYP2C9 African-American Coding;
		if Race=1 and z<=&AAprob6homo. then
			x1=1;

		if Race=1 and x1=0 and z1<=&AAprob3homo. then
			x2=1;

		if Race=1 and x1=0 and x2=0 and z2<=&AAprob5homo. then
			x3=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and z3<=&AAprob3and6. then
			x4=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and z4<=&AAprob11homo. then
			x5=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and 
			z5<=&AAprob5and6. then
				x6=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and 
			z6<=&AAprob11and6. then
				x7=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			z7<=&AAprob3and5. then
				x8=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and z8<=&AAprob11and3. then
				x9=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and z9<=&AAprob2and6. then
				x10=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and z10<=&AAprob11and5. then
				x11=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and z11<=&AAprob2homo. then
				x12=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and z12<=&AAprob2and3. then
				x13=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and 
			z13<=&AAprob2and5. then
				x14=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and 
			z14<=&AAprob11and2. then
				x15=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and z15<=&AAprob6and8. then
				x16=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and z16<=&AAprob3and8. then
				x17=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and z17<=&AAprob5and8. then
				x18=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and z18<=&AAprob11and8. then
				x19=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and z19<=&AAprob2and8. then
				x20=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and z20<=&AAprob8homo. then
				x21=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and x21=0 and 
			z21<=&AAprob1and6. then
				x22=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and x21=0 and x22=0 and 
			z22<=&AAprob1and3. then
				x23=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and x21=0 and x22=0 and 
			x23=0 and z23<=&AAprob1and5. then
				x24=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and x21=0 and x22=0 and 
			x23=0 and x24=0 and z24<=&AAprob1and11. then
				x25=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and x21=0 and x22=0 and 
			x23=0 and x24=0 and x25=0 and z25<=&AAprob1and2. then
				x26=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and x21=0 and x22=0 and 
			x23=0 and x24=0 and x25=0 and x26=0 and z26<=&AAprob1and8. then
				x27=1;

		if Race=1 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and x21=0 and x22=0 and 
			x23=0 and x24=0 and x25=0 and x26=0 and x27=0 then
				x28=1;

*CYP2C9 European Coding;
		if Race=0 and z<=&EUROprob6homo. then
			x1=1;

		if Race=0 and x1=0 and z1<=&EUROprob3homo. then
			x2=1;

		if Race=0 and x1=0 and x2=0 and z2<=&EUROprob5homo. then
			x3=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and z3<=&EUROprob3and6. then
			x4=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and z4<=&EUROprob11homo. then
			x5=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and 
			z5<=&EUROprob5and6. then
				x6=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and 
			z6<=&EUROprob11and6. then
				x7=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			z7<=&EUROprob3and5. then
				x8=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and z8<=&EUROprob11and3. then
				x9=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and z9<=&EUROprob2and6. then
				x10=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and z10<=&EUROprob11and5. then
				x11=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and z11<=&EUROprob2homo. then
				x12=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and z12<=&EUROprob2and3. then
				x13=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and 
			z13<=&EUROprob2and5. then
				x14=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and 
			z14<=&EUROprob11and2. then
				x15=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and z15<=&EUROprob6and8. then
				x16=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and z16<=&EUROprob3and8. then
				x17=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and z17<=&EUROprob5and8. then
				x18=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and z18<=&EUROprob11and8. then
				x19=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and z19<=&EUROprob2and8. then
				x20=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and z20<=&EUROprob8homo. then
				x21=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and x21=0 and 
			z21<=&EUROprob1and6. then
				x22=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and x21=0 and x22=0 and 
			z22<=&EUROprob1and3. then
				x23=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and x21=0 and x22=0 and 
			x23=0 and z23<=&EUROprob1and5. then
				x24=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and x21=0 and x22=0 and 
			x23=0 and x24=0 and z24<=&EUROprob1and11. then
				x25=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and x21=0 and x22=0 and 
			x23=0 and x24=0 and x25=0 and z25<=&EUROprob1and2. then
				x26=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and x21=0 and x22=0 and 
			x23=0 and x24=0 and x25=0 and x26=0 and z26<=&EUROprob1and8. then
				x27=1;

		if Race=0 and x1=0 and x2=0 and x3=0 and x4=0 and x5=0 and x6=0 and x7=0 and 
			x8=0 and x9=0 and x10=0 and x11=0 and x12=0 and x13=0 and x14=0 and x15=0 
			and x16=0 and x17=0 and x18=0 and x19=0 and x20=0 and x21=0 and x22=0 and 
			x23=0 and x24=0 and x25=0 and x26=0 and x27=0 then
				x28=1;
		*VKORC1 A/G(x5);
		probAAMinorAlleleHomo=&AAVKORC1.;
		probAAMajorAlleleHomo=1-&AAVKORC1.;
		probEUMinorAlleleHomo=&EUROVKORC1.;
		probEUMajorAlleleHomo=1-&EUROVKORC1.;

		z30=ranuni(&i+739281079217);

		VKORC=1;

		if Race = 1 and z30 >=1-probAAMajorAlleleHomo then
			VKORC=0;
		if Race = 0 and z30 >=1-probEUMajorAlleleHomo then 
			VKORC=0;

		if Race = 1 and z30 <=probAAMinorAlleleHomo then
			VKORC=2;
		if Race = 0 and z30 <=probEUMinorAlleleHomo then 
			VKORC=2;

		if VKORC=1 then
			VKORCAG=1;
		else
			VKORCAG=0;

		if VKORC=2 then
			VKORCAA=1;
		else
			VKORCAA=0;

		*Amiodarone status of population(Amiodarone);

		if ranuni(&i+137587) le &percentAMIO. then
			Amiodarone=1;
		else
			Amiodarone=0;
		*1=taking amidarone and 0=not taking amiodarone
		*dosage formula (no changes);

		*Dosing Formula if Genes Unknown;
		nogeneweeklydosage=((5.6044-
		0.2614*Age+
		0.0087*Height+
		0.0128*Weight-
		0.5503*Amiodarone-
		0.4854-
		0.2188-
		0.2760*Race)**2);

		*Dosing Formula if Genes Known;
		geneweeklydosage=((5.6044-
		0.2614*Age+
		0.0087*Height+
		0.0128*Weight-
		0.5503*Amiodarone-
		0.8677*VKORCAG-
		1.6974*VKORCAA-
		0.5211*x26-
		0.9357*x23-
		1.0616*x12-
		1.9206*x13-
		2.3312*x2-
		0.2760*Race)**2);

		*Dosing Formula if AA Genes Known;
		AAgeneweeklydosage = geneweeklydosage;
			
			if x1=1 or x3=1 or x5=1 or x6=1 or x7=1 or x11=1 or x16=1 or x18=1 or x19=1 or x21=1 then 
				AAgeneweeklydosage = geneweeklydosage*0.7;
			if x4=1 or x8=1 or x9=1 or x10=1 or x14=1 or x15=1 or x17=1 or x20=1 or x22=1 or x24=1 or x25=1 or x27=1 then 
				AAgeneweeklydosage = geneweeklydosage*0.85;

	*Overdosed Patients;
		If nogeneweeklydosage <= 21 then Clinical = "Low";
		If 21 < nogeneweeklydosage < 49 then Clinical = "Med";
		If nogeneweeklydosage >= 49 then Clinical = "Hi";

		If geneweeklydosage <= 21 then Pgx = "Low";
		If 21 < geneweeklydosage < 49 then Pgx = "Med";
		If geneweeklydosage >= 49 then Pgx = "Hi";

		If AAgeneweeklydosage <= 21 then AAPgx = "Low";
		If 21 < AAgeneweeklydosage < 49 then AAPgx = "Med";
		If AAgeneweeklydosage >= 49 then AAPgx = "Hi";

		z31 = ranuni(&i+40000000);
		z32 = ranuni(&i+50000000);
		z33 = ranuni(&i+60000000);

		If x1=1 or x3=1 or x5=1 or x6=1 or x7=1 or x11=1 or x16=1 or x18=1 or x19=1 or x21=1 or x4=1 or x8=1 
			or x9=1 or x10=1 or x14=1 or x15=1 or x17=1 or x20=1 or x22=1 or x24=1 or x25=1 or x27=1 
			then True = AAPgx;

			else if pgx = "Hi" and z31 le &probhightruepgx. then True = Pgx;
				else if pgx = "Hi" and z31 > &probhightruepgx. then True = "Med";
			else if pgx = "Med" and z32 le &probmedoverdosedpgx. then True = "Low";
				else if pgx = "Med" and &probmedoverdosedpgx. < z32 <= &probmedtruepgx.
					then True = Pgx;
				else if pgx = "Med" and z32 < &probmedtruepgx. then True = "Hi";
			else if pgx = "Low" and z33 le &problowtruepgx. then True = Pgx;
				else if pgx = "Low" and z33 > &problowtruepgx. then True = "Med";



		ClinicalOverdosed =0;
			If true = "Med" and Clinical = "Hi" then ClinicalOverdosed = 1;
			If true = "Low" and Clinical = "Hi" then ClinicalOverdosed = 1;
			If true = "Low" and Clinical = "Med" then ClinicalOverdosed = 1;
		GeneticsOverdosed =0;
			If true = "Med" and Pgx = "Hi" then GeneticsOverdosed = 1;
			If true = "Low" and Pgx = "Hi" then GeneticsOverdosed = 1;
			If true = "Low" and Pgx = "Med" then GeneticsOverdosed = 1;
		AAGeneticsOverdosed =0;
			If true = "Med" and AAPgx = "Hi" then AAGeneticsOverdosed = 1;
			If true = "Low" and AAPgx = "Hi" then AAGeneticsOverdosed = 1;
			If true = "Low" and AAPgx = "Med" then AAGeneticsOverdosed = 1;


		ClinicalUnderdosed =0;
			If true = "Hi" and Clinical = "Med" then ClinicalUnderdosed = 1;
			If true = "Hi" and Clinical = "Low" then ClinicalUnderdosed = 1;
			If true = "Med" and Clinical = "Hi" then ClinicalUnderdosed = 1;
		GeneticsUnderdosed =0;
			If true = "Hi" and Pgx = "Med" then GeneticsUnderdosed = 1;
			If true = "Hi" and Pgx = "Low" then GeneticsUnderdosed = 1;
			If true = "Med" and Pgx = "Hi" then GeneticsUnderdosed = 1;
		AAGeneticsUnderdosed =0;
			If true = "Hi" and AAPgx = "Med" then AAGeneticsUnderdosed = 1;
			If true = "Hi" and AAPgx = "Low" then AAGeneticsUnderdosed = 1;
			If true = "Med" and AAPgx = "Hi" then AAGeneticsUnderdosed = 1;

		output temp1;
	end;
run;

ods  output Summary=MeansStatistics;

proc means data=temp1 min max mean std range median q1 q3;
	class race;
	var age weight height nogeneweeklydosage geneweeklydosage AAgeneweeklydosage;
	run;


proc sort data=temp1;
	by race;
	run;

proc freq data=temp1;
	by race;
	table x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 
		x24 x25 x26 x27 x28 VKORC Amiodarone clinical pgx aapgx true clinicaloverdosed geneticsoverdosed
		aageneticsoverdosed clinicalunderdosed geneticsunderdosed aageneticsunderdosed;
	run;

