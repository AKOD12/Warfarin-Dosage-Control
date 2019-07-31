options nonotes;
options nodate formdlim="-";

%macro simulation(scenar,minorallelefreq,drugadr,qtdrug,adrqt,qtsnp,qtdrugsnp);

/********Static Variables*****************/
%Let nsim=100000;       		*number of simulatiosn;
%Let nstudy=120000;    		*study population size;
%Let meanage=54.4;			*mean age;
%Let sdage=5.7;				*standard deviation of age;
%Let percentmale=0.473;		*percent male;
%Let meanU=0;				*mean unmeasured/unknown covariates;
%Let stU=4.95;				*standard deviation of U;
%Let artintercept=-6.00;	*intercept in art derivation;
%Let artbage=0.075;			*alpha of age in art derivation;
%Let htnintercept=-5.25;	*intercept in HTN1 derivation;
%Let htnage=0.08;			*alpha of age in HTN1 derivation;
%Let htnsex=0.16;			*alpha of sex in HTN1 derivation;
%Let htn2intercept=-4.74;	*intercept in HTN2-4 derivation;
%Let htn2age=0.05;			*alpha of age in HTN2-4 derivation;
%Let htn2sex=0.19;			*alpha of sex in HTN2-4 derivation;
%Let htn2htn=3.92;			*alpha of HTNi-1 in HTN2-4 derivation;
%Let probtreatment=0.91;	*probabilty of HTN treatment;
%Let drugintercept=1.55;	*intercept in drug derivation;
%Let drugart=-90;			*alpha of art in drug derivation;
%Let drugdrug=4.80;			*alpha of drugi-1 in drug derivation;
/**********%Let drugadr=-50; *-3.75;			*alpha of adr in drug derivation;**/
%Let qtintercept=389;		*mean QT;
/*********%Let qtdrug=30;  *5;				*beta of drug in qt derivation;********/
%Let qtage=0.39;			*beta of age in qt derivation;
%Let qtsex=-3.84;			*beta of sex in qt derivation;
%Let qtU=4.35;				*beta of U in qt derivation;
%Let qtart=-4.90;			*beta of art in qt derivation;
%Let qthtn=4.50;			*beta of HTN in qt derivation;
%Let stqt=14;				*standard error of qt derivation;
%Let adrintercept=-4.55;	*intercept of ADR derivation;
%Let adrsnp=0;			*alpha of SNP in ADR derivation;
/********%Let adrqt=50;*2.05;			*alpha of QTlong in ADR derivation;********/
%Let lossintercept=-2.62;	*intercept in loss-to-follow-up derivation;
%Let lossage=-0.01;			*alpha of age in loss-to-follow-up derivation;
%Let lossadr=2.89;			*aloph of ADR in loss-to-follow-up derivation;

/*****************************************/
/******DATA GENERATION********************/
/*****************************************/

%do i=1 %to &nsim.;  
	data temp1;
	do j=1 to &nstudy.;

	*Coding Genotype (x1);
	probMinorAlleleHomo = &minorallelefreq.**2;
	probMajorAlleleHomo = (1-&minorallelefreq.)**2;

	z=ranuni(&i+2000000);

	x1 = 1;
	if z >= 1-probMajorAlleleHomo then x1 = 0;
	if z <= probMinorAlleleHomo then x1 = 2;

	*Coding Age (x2);
	x2 = &meanage. + &sdage. * rannor(&i+3000000);

	*Coding Gender (x3);
	if ranuni(&i+4000000) le &percentmale. then x3 = 1;
	else x3 = 0;

	*Coding U (x4);
	x4 = &meanU. + &stU. * rannor(&i+5000000);

	*Coding Diabetes (x5);
	pDiab = (1 + exp(-(&diabintercept. + &diabage.*x2)))**-1;
	if ranuni(&i+7000000) le pDiab then x5 = 1;
	else x5 = 0;

	*Coding HTN1 (x6_1);
	pHTN1 = (1 + exp(-(&htnintercept. + &htnage.*x2 + &htnsex.*x3)))**-1;
	if ranuni(&i+8000000) le pHTN1 then x6_1 = 1;
	else x6_1 = 0;

	*Coding HTN2 (x6_2);
	pHTN2 = (1 + exp(-(&htn2intercept. + &htn2age.*x2 + &htn2sex.*x3 + &htn2htn.*x6_1)))**-1;
	if ranuni(&i+9000000) le pHTN2 then x6_2 = 1;
	else x6_2 = 0;

	*Coding HTN3 (x6_3);
	pHTN3 = (1 + exp(-(&htn2intercept. + &htn2age.*x2 + &htn2sex.*x3 + &htn2htn.*x6_2)))**-1;
	if ranuni(&i+10000000) le pHTN3 then x6_3 = 1;
	else x6_3 = 0;

	*Coding HTN4 (x6_4);
	pHTN4 = (1 + exp(-(&htn2intercept. + &htn2age.*x2 + &htn2sex.*x3 + &htn2htn.*x6_3)))**-1;
	if ranuni(&i+11000000) le pHTN4 then x6_4 = 1;
	else x6_4 = 0;

	*Coding Treatment1 (x7_1);
	/**if x6_1 = 0 then pTreat = 0;
	else if x6_1 = 1 then pTreat = (1 + exp(-(&probtreatment. + rannor(&i+9000000))))**-1;
	if ranuni(&i+12000000) le pTreat then x7_1 = 1;
	else x7_1 = 0;**/

	*Coding Treatment2 (x7_2);
	if x6_2 = 0 then pTreat2 = 0;
	else if x6_2 = 1 then pTreat2 = (1 + exp(-(&probtreatment. + rannor(&i+9000000))))**-1;
	if ranuni(&i+13000000) le pTreat2 then x7_2 = 1;
	else x7_2 = 0;

	*Coding Treatment3 (x7_3);
	if x6_3 = 0 then pTreat3 = 0;
	if x6_3 = 1 then pTreat3 = (1 + exp(-(&probtreatment. + rannor(&i+9000000))))**-1;
	if ranuni(&i+14000000) le pTreat3 then x7_3 = 1;
	else x7_3 = 0;

	*Coding Treatment4 (x7_4);
	if x6_4 = 0 then pTreat4 = 0;
	if x6_4 = 1 then pTreat4 = (1 + exp(-(&probtreatment. + rannor(&i+9000000))))**-1;
	if ranuni(&i+15000000) le pTreat4 then x7_4 = 1;
	else x7_4 = 0;

	*Coding Drug1 (x8_1) and Comparator1 (X9_1);
	x8_1 = 0;
	x9_1 = .;

	*Coding QT1 (QT_1)  ;
	ExpQT_1 = (&QTintercept. + &qtsnp.*x1 + &qtdrug.*x8_1 + &qtdrugsnp.*x1*x8_1 + &qtage.*x2 + &qtsex.*x3 + &qtU.*x4 + &qtdiab.*x5 + &qthtn.*x6_1);
	QT_1 = ExpQT_1 + &stqt.*rannor(&i+18000000);
	if QT_1 >= 450 then QTL_1 = 1;
	else QTL_1 = 0;

	*Coding Drug2 (x8_2) and Comparator2 (X9_2) + &drugadr.*x10_2;
	if x7_2 = 0 then pDrug2 = 0;
	else if x7_2 = 1 then pDrug2 = (1 + exp(-(&drugintercept. + &drugdiab.*x5 + &drugdrug.*x8_1)))**-1;
	if x7_2 = 0 then pComparator2 = 0;
	else pComparator2 = 1 - pDrug2;
	if ranuni(&i+20000000) le pDrug2 then x8_2 = 1;
	else x8_2 = 0;
	if ranuni(&i+21000000) le pComparator2 then x9_2 = 1;
	else x9_2 = 0;

	*Coding QT2 (QT_2);
	ExpQT_2 = (&QTintercept. + &qtdrug.*x8_2 + &qtsnp.*x1 + &qtdrugsnp.*x1*x8_2 + &qtage.*x2 + &qtsex.*x3 + &qtU.*x4 + &qtdiab.*x5 + &qthtn.*x6_2);
	QT_2 = ExpQT_2 + &stqt.*rannor(&i+22000000);
	if QT_2 >= 450 then QTL_2 = 1;
	else QTL_2 = 0;

	*Coding ADR3 (X10_3);
	if x8_2 = 0 then pADR3 = 0;
	else pADR3 = (1 + exp(-(&adrintercept. + &adrsnp.*x1 + &adrqt.*qtl_2)))**-1;
	if ranuni(&i+23000000) le pADR3 then x10_3 = 1;
	else x10_3 = 0; 

	*Coding Drug3 (x8_3) and Comparator2 (X9_3);
	if x7_3 = 0 then pDrug3 = 0;
	else if x7_3 = 1 then pDrug3 = (1 + exp(-(&drugintercept. + &drugdiab.*x5 + &drugdrug.*x8_2 + &drugadr.*x10_3)))**-1;
	if x7_3 = 0 then pComparator3 = 0;
	else pComparator3 = 1 - pDrug3;
	if ranuni(&i+24000000) le pDrug3 then x8_3 = 1;
	else x8_3 = 0;
	if ranuni(&i+25000000) le pComparator3 then x9_3 = 1;
	else x9_3 = 0;

	*Coding QT3 (QT_3);
	ExpQT_3 = (&QTintercept. + &qtdrug.*x8_3 + &qtsnp.*x1 + &qtdrugsnp.*x1*x8_3 + &qtage.*x2 + &qtsex.*x3 + &qtU.*x4 + &qtdiab.*x5 + &qthtn.*x6_3);
	QT_3 = ExpQT_3 + &stqt.*rannor(&i+26000000);
	if QT_3 >= 450 then QTL_3 = 1;
	else QTL_3 = 0;

	*Coding ADR4 (X10_4);
	if x8_3 = 0 then pADR4 = 0;
	else pADR4 = (1 + exp(-(&adrintercept. + &adrsnp.*x1 + &adrqt.*qtl_3)))**-1;
	if ranuni(&i+27000000) le pADR4 then x10_4 = 1;
	else x10_4 = 0; 

	*Coding Drug4 (x8_4) and Comparator2 (X9_4);
	if x7_4 = 0 then pDrug4 = 0;
	else if x7_4 = 1 then pDrug4 = (1 + exp(-(&drugintercept. + &drugdiab.*x5 + &drugdrug.*x8_3 + &drugadr.*x10_4)))**-1;
	if x7_4 = 0 then pComparator4 = 0;
	else pComparator4 = 1 - pDrug4;
	if ranuni(&i+28000000) le pDrug4 then x8_4 = 1;
	else x8_4 = 0;
	if ranuni(&i+29000000) le pComparator4 then x9_4 = 1;
	else x9_4 = 0;

	*Coding QT4 (QT_4);
	ExpQT_4 = (&QTintercept. + &qtdrug.*x8_4 + &qtsnp.*x1 + &qtdrugsnp.*x1*x8_4 + &qtage.*x2 + &qtsex.*x3 + &qtU.*x4 + &qtdiab.*x5 + &qthtn.*x6_4);
	QT_4 = ExpQT_4 + &stqt.*rannor(&i+30000000);

	*Coding Loss-to-Follow-Up3 (L_3);
	pLoss3 = (1 + exp(-(&lossintercept. + &lossage.*x2 + &lossadr.*x10_3)))**-1;
	if ranuni(&i+32000000) le pLoss3 then L_3 = 1;
	else L_3 = 0;

	*Coding Loss-to-Follow-Up2 (L_2);
	pLoss4 = (1 + exp(-(&lossintercept. + &lossage.*x2 + &lossadr.*x10_4)))**-1;
	if ranuni(&i+33000000) le pLoss4 then L_4 = 1;
	else L_4 = 0;

	scenario=&scenar.;

	one = 1;

	x10_1 = .;
	L_1 = .;

	output temp1;

	end;
	run;

	*Applying Loss-to-Followup;
	data temp1_loss;
	set temp1;
	if L_3 = 1 then do;
		x6_3 = .; x6_4 = .;
		x7_3 = .; x7_4 = .;
		x8_3 = .; x8_4 = .;
		x9_3 = .; x9_4 = .;
		x10_3 = .; x10_4 = .;
		qt_3 = .; qt_4 = .;
		end;
	if L_4 = 1 then do;
		x6_4 = .;
		x7_4 = .;
		x8_4 = .;
		x9_4 = .;
		x10_4 = .;
		qt_4 = .;
		end;
	x10_2 = .;
	l_2 = .;
	run;

	*Converting Data to Long Format For Analysis;
	data temp1_long;
	set temp1_loss;
	array htn (4) x6_1 x6_2 x6_3 x6_4;
	array drug (4) x8_1 x8_2 x8_3 x8_4;
	array comp (4) x9_1 x9_2 x9_3 x9_4;
	array qtint (4) qt_1 qt_2 qt_3 qt_4;
	array adr (4) x10_1 x10_2 x10_3 x10_4;
	array loss (4) L_1 L_2 L_3 L_4;
	do inventory = 1 to 4;
		x6 = htn(inventory);
		x8 = drug(inventory);
		x9 = comp(inventory);
		qt = qtint(inventory);
		x10 = adr(inventory);
		L = loss(inventory);
		output;
	end;
	run;


/*****************************************/
/******MODEL SPECIFICATION****************/
/*****************************************/

/******UNBIASED RESULT********************/

ods noresults;
ods trace on;
ods output ParameterEstimates=res;

	proc genmod data=temp1 desc;
	model qt_2 = x1 x8_2 x8_2*x1 x2 x3 x5 x6_2 / dist=normal link=identity;
	title "Cross-Sectional New User Unbiased";
	run;

	data unbiased_a;
	set res;
	if parameter = 'x8_2';
	ueit = estimate;
	vueit = stderr**2;
	if probchisq < 0.00000005 then sueit = 1;
		else sueit = 0;
	if lowerwaldcl < &qtdrug. < upperwaldcl then covueit = 1;
		else covueit = 0;
	biasueit = ueit - &qtdrug.;
	squeit = (ueit - &qtdrug.)**2;
	msqueit = vueit + squeit;
	ueprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;
	
	data unbiased_b;
	set res;
	if parameter = 'x1*x8_2';
	uintit = estimate;
	vuintit = stderr**2;
	if probchisq < 0.00000005 then suintit = 1; 
		else suintit = 0;
	if lowerwaldcl < &qtdrugsnp. < upperwaldcl then covuintit = 1;
		else covuintit = 0;
	biasuintit = uintit - &qtdrugsnp.;
	squintit = (uintit - &qtdrugsnp.)**2;
	msquintit = vuintit + squintit;
	uintprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;

	data unbiased_c;
	set res;
	if parameter = 'x1';
	ugit = estimate;
	vugit = stderr**2;
	if probchisq < 0.00000005 then sugit = 1; 
		else sugit = 0;
	if lowerwaldcl < &qtsnp. < upperwaldcl then covugit = 1;
		else covugit = 0;
	biasugit = ugit - &qtsnp.;
	squgit = (ugit - &qtsnp.)**2;
	msqugit = vugit + squgit;
	ugprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;


 /******UNBIASED RESULT: ACTIVE COMPARATOR*****/

ods noresults;
ods trace on;
ods output ParameterEstimates=res;

	data temp1_nuac;
	set temp1;
	if x8_2 = 1 then keep = 1;
	if x9_2 = 1 then keep = 1;
		else keep = 0;
	if keep = 1;
	drop keep;
	run;

	proc genmod data=temp1_nuac desc;
	model qt_2 = x1 x8_2 x8_2*x1 x2 x3 x5 x6_2 / dist=normal link=identity;
	title "Cross-Sectional New User Unbiased: Active Comparator";
	run;

	data unbiased_ac_a;
	set res;
	if parameter = 'x8_2';
	uaceit = estimate;
	vuaceit = stderr**2;
	if probchisq < 0.00000005 then suaceit = 1;
		else suaceit = 0;
	if lowerwaldcl < &qtdrug. < upperwaldcl then covuaceit = 1;
		else covuaceit = 0;
	biasuaceit = uaceit - &qtdrug.;
	squaceit = (uaceit - &qtdrug.)**2;
	msquaceit = vuaceit + squaceit;
	uaceprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;
	
	data unbiased_ac_b;
	set res;
	if parameter = 'x1*x8_2';
	uacintit = estimate;
	vuacintit = stderr**2;
	if probchisq < 0.00000005 then suacintit = 1; 
		else suacintit = 0;
	if lowerwaldcl < &qtdrugsnp. < upperwaldcl then covuacintit = 1;
		else covuacintit = 0;
	biasuacintit = uacintit - &qtdrugsnp.;
	squacintit = (uacintit - &qtdrugsnp.)**2;
	msquacintit = vuacintit + squacintit;
	uacintprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;

	data unbiased_ac_c;
	set res;
	if parameter = 'x1';
	uacgit = estimate;
	vuacgit = stderr**2;
	if probchisq < 0.00000005 then suacgit = 1; 
		else suacgit = 0;
	if lowerwaldcl < &qtsnp. < upperwaldcl then covuacgit = 1;
		else covuacgit = 0;
	biasuacgit = uacgit - &qtsnp.;
	squacgit = (uacgit - &qtsnp.)**2;
	msquacgit = vuacgit + squacgit;
	uacgprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;

/******LONGITUDINAL ANALYSES**************/

/******WHOLE COHORT***********************/
ods noresults;
ods trace on;
ods output GEEEmpPEst=res;

	data temp1_longwc;
	set temp1_long;
	if inventory = 1 then delete;
	if inventory = 2 then delete;
	run;

	proc genmod data=temp1_longwc desc;
	class j;
	model qt = x1 x8 x8*x1 x2 x3 x5 x6 / dist=normal link=identity;
	repeated subject = j / type=ind;
	title "Longitudinal Data: Whole Cohort";
	run; quit;

	data long_wc_a;
	set res;
	if parm = 'x8';
	eit = estimate;
	veit = stderr**2;
	if probz < 0.00000005 then seit = 1;
		else seit = 0;
	if lowercl < &qtdrug. < uppercl then coveit = 1;
		else coveit = 0;
	biaseit = eit - &qtdrug.;
	sqeit = (eit - &qtdrug.)**2;
	msqeit = veit + sqeit;
	eprobz = probz;
	drop parm estimate stderr lowercl uppercl probz z;
	run;
	
	data long_wc_b;
	set res;
	if parm = 'x1*x8';
	intit = estimate;
	vintit = stderr**2;
	if probz < 0.00000005 then sintit = 1; 
		else sintit = 0;
	if lowercl < &qtdrugsnp. < uppercl then covintit = 1;
		else covintit = 0;
	biasintit = intit - &qtdrugsnp.;
	sqintit = (intit - &qtdrugsnp.)**2;
	msqintit = vintit + sqintit;
	intprobz = probz;
	drop parm estimate stderr lowercl uppercl probz z;
	run;

	data long_wc_c;
	set res;
	if parm = 'x1';
	git = estimate;
	vgit = stderr**2;
	if probz < 0.00000005 then sgit = 1; 
		else sgit = 0;
	if lowercl < &qtsnp. < uppercl then covgit = 1;
		else covgit = 0;
	biasgit = git - &qtsnp.;
	sqgit = (git - &qtsnp.)**2;
	msqgit = vgit + sqgit;
	gprobz = probz;
	drop parm estimate stderr lowercl uppercl probz z;
	run;


/******ACTIVE COMPARATOR******************/
ods noresults;
ods output GEEEmpPEst=res;

	data temp1_ac;
	set temp1_longwc;
	if x8 = 1 then keep = 1;
	if x9 = 1 then keep = 1;
		else keep = 0;
	if keep = 1;
	drop keep;
	run;

	proc genmod data=temp1_ac desc;
	class j;
	model qt = x1 x8 x8*x1 x2 x3 x5 x6 / dist=normal link=identity;
	repeated subject = j / type=ind;
	title "Longitudinal Data: Active Comparator";
	run; quit;

	data long_ac_a;
	set res;
	if parm = 'x8';
	aceit = estimate;
	vaceit = stderr**2;
	if probz < 0.05 then saceit = 1;
		else saceit = 0;
	if lowercl < &qtdrug. < uppercl then covaceit = 1;
		else covaceit = 0;
	biasaceit = aceit - &qtdrug.;
	sqaceit = (aceit - &qtdrug.)**2;
	msqaceit = vaceit + sqaceit;
	aceprobz = probz;
	drop parm estimate stderr lowercl uppercl probz z;
	run;
	
	data long_ac_b;
	set res;
	if parm = 'x1*x8';
	acintit = estimate;
	vacintit = stderr**2;
	if probz < 0.05 then sacintit = 1; 
		else sacintit = 0;
	if lowercl < &qtdrugsnp. < uppercl then covacintit = 1;
		else covacintit = 0;
	biasacintit = acintit - &qtdrugsnp.;
	sqacintit = (acintit - &qtdrugsnp.)**2;
	msqacintit = vacintit + sqacintit;
	acintprobz = probz;
	drop parm estimate stderr lowercl uppercl probz z;
	run;

	data long_ac_c;
	set res;
	if parm = 'x1';
	acgit = estimate;
	vacgit = stderr**2;
	if probz < 0.05 then sacgit = 1; 
		else sacgit = 0;
	if lowercl < &qtsnp. < uppercl then covacgit = 1;
		else covacgit = 0;
	biasacgit = acgit - &qtsnp.;
	sqacgit = (acgit - &qtsnp.)**2;
	msqacgit = vacgit + sqacgit;
	acgprobz = probz;
	drop parm estimate stderr lowercl uppercl probz z;
	run;


/******CROSS-SECTIONAL************/

	data temp1_cs;
	set temp1_longwc;
	if inventory = 4 then keep = 1;
		else keep = 0;
	if keep = 1;
	run;


/******WHOLE COHORT***********************/
ods noresults;
ods output ParameterEstimates=res;

	proc genmod data=temp1_cs desc;
	model qt = x1 x8 x8*x1 x2 x3 x5 x6 / dist=normal link=identity;
	title "Cross-Sectional: Whole Cohort";
	run; quit;

	data cs_wc_a;
	set res;
	if parameter = 'x8';
	cseit = estimate;
	vcseit = stderr**2;
	if probchisq < 0.00000005 then scseit = 1;
		else scseit = 0;
	if lowerwaldcl < &qtdrug. < upperwaldcl then covcseit = 1;
		else covcseit = 0;
	biascseit = cseit - &qtdrug.;
	sqcseit = (cseit - &qtdrug.)**2;
	msqcseit = vcseit + sqcseit;
	cseprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;

	data cs_wc_b;
	set res;
	if parameter = 'x1*x8';
	csintit = estimate;
	vcsintit = stderr**2;
	if probchisq < 0.00000005 then scsintit = 1;
		else scsintit = 0;
	if lowerwaldcl < &qtdrugsnp. < upperwaldcl then covcsintit = 1;
		else covcsintit = 0;
	biascsintit = csintit - &qtdrugsnp.;
	sqcsintit = (csintit - &qtdrugsnp.)**2;
	msqcsintit = vcsintit + sqcsintit;
	csintprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;
	
	data cs_wc_c;
	set res;
	if parameter = 'x1';
	csgit = estimate;
	vcsgit = stderr**2;
	if probchisq < 0.00000005 then scsgit = 1;
		else scsgit = 0;
	if lowerwaldcl < &qtsnp. < upperwaldcl then covcsgit = 1;
		else covcsgit = 0;
	biascsgit = csgit - &qtsnp.;
	sqcsgit = (csgit - &qtsnp.)**2;
	msqcsgit = vcsgit + sqcsgit;
	csgprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;


/******ACTIVE COMPARATOR******************/
ods noresults;
ods output ParameterEstimates=res;

	data temp1_ac_cs;
	set temp1_cs;
	if x8 = 1 then keep = 1;
	if x9 = 1 then keep = 1;
		else keep = 0;
	if keep = 1;
	drop keep;
	run;

	proc genmod data=temp1_ac_cs desc;
	model qt = x1 x8 x8*x1 x2 x3 x5 x6 / dist=normal link=identity;
	title "Cross-Sectional: Active Comparator";
	run; quit;

	data cs_ac_a;
	set res;
	if parameter = 'x8';
	csaceit = estimate;
	vcsaceit = stderr**2;
	if probchisq < 0.00000005 then scsaceit = 1;
		else scsaceit = 0;
	if lowerwaldcl < &qtdrug. < upperwaldcl then covcsaceit = 1;
		else covcsaceit = 0;
	biascsaceit = csaceit - &qtdrug.;
	sqcsaceit = (csaceit - &qtdrug.)**2;
	msqcsaceit = vcsaceit + sqcsaceit;
	csaceprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;
	
	data cs_ac_b;
	set res;
	if parameter = 'x1*x8';
	csacintit = estimate;
	vcsacintit = stderr**2;
	if probchisq < 0.00000005 then scsacintit = 1;
		else scsacintit = 0;
	if lowerwaldcl < &qtdrugsnp. < upperwaldcl then covcsacintit = 1;
		else covcsacintit = 0;
	biascsacintit = csacintit - &qtdrugsnp.;
	sqcsacintit = (csacintit - &qtdrugsnp.)**2;
	msqcsacintit = vcsacintit + sqcsacintit;
	csacintprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;
	
	data cs_ac_c;
	set res;
	if parameter = 'x1';
	csacgit = estimate;
	vcsacgit = stderr**2;
	if probchisq < 0.00000005 then scsacgit = 1;
		else scsacgit = 0;
	if lowerwaldcl < &qtsnp. < upperwaldcl then covcsacgit = 1;
		else covcsacgit = 0;
	biascsacgit = csacgit - &qtsnp.;
	sqcsacgit = (csacgit - &qtsnp.)**2;
	msqcsacgit = vcsacgit + sqcsacgit;
	csacgprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;

/******CROSS-SECTIONAL VERSION 2************/

	data temp1_cs2;
	set temp1_longwc;
	if inventory = 3 then keep = 1;
		else keep = 0;
	if keep = 1;
	run;


/******WHOLE COHORT***********************/
ods noresults;
ods output ParameterEstimates=res;

	proc genmod data=temp1_cs2 desc;
	model qt = x1 x8 x8*x1 x2 x3 x5 x6 / dist=normal link=identity;
	title "Cross-Sectional Version 2: Whole Cohort";
	run; quit;

	data cs2_wc_a;
	set res;
	if parameter = 'x8';
	cs2eit = estimate;
	vcs2eit = stderr**2;
	if probchisq < 0.00000005 then scs2eit = 1;
		else scs2eit = 0;
	if lowerwaldcl < &qtdrug. < upperwaldcl then covcs2eit = 1;
		else covcs2eit = 0;
	biascs2eit = cs2eit - &qtdrug.;
	sqcs2eit = (cs2eit - &qtdrug.)**2;
	msqcs2eit = vcs2eit + sqcs2eit;
	cs2eprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;

	data cs2_wc_b;
	set res;
	if parameter = 'x1*x8';
	cs2intit = estimate;
	vcs2intit = stderr**2;
	if probchisq < 0.00000005 then scs2intit = 1;
		else scs2intit = 0;
	if lowerwaldcl < &qtdrugsnp. < upperwaldcl then covcs2intit = 1;
		else covcs2intit = 0;
	biascs2intit = cs2intit - &qtdrugsnp.;
	sqcs2intit = (cs2intit - &qtdrugsnp.)**2;
	msqcs2intit = vcs2intit + sqcs2intit;
	cs2intprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;
	
	data cs2_wc_c;
	set res;
	if parameter = 'x1';
	cs2git = estimate;
	vcs2git = stderr**2;
	if probchisq < 0.00000005 then scs2git = 1;
		else scs2git = 0;
	if lowerwaldcl < &qtsnp. < upperwaldcl then covcs2git = 1;
		else covcs2git = 0;
	biascs2git = cs2git - &qtsnp.;
	sqcs2git = (cs2git - &qtsnp.)**2;
	msqcs2git = vcs2git + sqcs2git;
	cs2gprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;


/******ACTIVE COMPARATOR******************/
ods noresults;
ods output ParameterEstimates=res;

	data temp1_ac_cs2;
	set temp1_cs2;
	if x8 = 1 then keep = 1;
	if x9 = 1 then keep = 1;
		else keep = 0;
	if keep = 1;
	drop keep;
	run;

	proc genmod data=temp1_ac_cs2 desc;
	model qt = x1 x8 x8*x1 x2 x3 x5 x6 / dist=normal link=identity;
	title "Cross-Sectional Version 2: Active Comparator";
	run; quit;

	data cs2_ac_a;
	set res;
	if parameter = 'x8';
	cs2aceit = estimate;
	vcs2aceit = stderr**2;
	if probchisq < 0.00000005 then scs2aceit = 1;
		else scs2aceit = 0;
	if lowerwaldcl < &qtdrug. < upperwaldcl then covcs2aceit = 1;
		else covcs2aceit = 0;
	biascs2aceit = cs2aceit - &qtdrug.;
	sqcs2aceit = (cs2aceit - &qtdrug.)**2;
	msqcs2aceit = vcs2aceit + sqcs2aceit;
	cs2aceprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;
	
	data cs2_ac_b;
	set res;
	if parameter = 'x1*x8';
	cs2acintit = estimate;
	vcs2acintit = stderr**2;
	if probchisq < 0.00000005 then scs2acintit = 1;
		else scs2acintit = 0;
	if lowerwaldcl < &qtdrugsnp. < upperwaldcl then covcs2acintit = 1;
		else covcs2acintit = 0;
	biascs2acintit = cs2acintit - &qtdrugsnp.;
	sqcs2acintit = (cs2acintit - &qtdrugsnp.)**2;
	msqcs2acintit = vcs2acintit + sqcs2acintit;
	cs2acintprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;
	
	data cs2_ac_c;
	set res;
	if parameter = 'x1';
	cs2acgit = estimate;
	vcs2acgit = stderr**2;
	if probchisq < 0.00000005 then scs2acgit = 1;
		else scs2acgit = 0;
	if lowerwaldcl < &qtsnp. < upperwaldcl then covcs2acgit = 1;
		else covcs2acgit = 0;
	biascs2acgit = cs2acgit - &qtsnp.;
	sqcs2acgit = (cs2acgit - &qtsnp.)**2;
	msqcs2acgit = vcs2acgit + sqcs2acgit;
	cs2acgprobz = probchisq;
	drop parameter estimate stderr lowerwaldcl upperwaldcl probchisq chisq;
	run;


/*****************************************/
/******COMBINE DATASETS*******************/
/*****************************************/

	data all;
	merge unbiased_a unbiased_b unbiased_c unbiased_ac_a unbiased_ac_b unbiased_ac_c 
			long_wc_a long_wc_b long_wc_c long_ac_a long_ac_b long_ac_c 
			cs_wc_a cs_wc_b cs_wc_c cs_ac_a cs_ac_b cs_ac_c
			cs2_wc_a cs2_wc_b cs2_wc_c cs2_ac_a cs2_ac_b cs2_ac_c;
	scenario = &scenar.;
	run;

	proc append base = results.completesp__&scenar.P1 data=all;
	run;

	DM log 'clear' wedit continue;
	DM output 'clear' wedit continue;
	run;

%end;

proc append base=results.datasets1__&scenar.P1 data=temp1;
run;



%mend simulation;

/******0.45*************/

%simulation(1,0.25,-50,30,50,0,0);

