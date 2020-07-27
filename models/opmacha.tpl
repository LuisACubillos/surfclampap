TOP_OF_MAIN_SECTION
  arrmblsize=300000; // 
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(30000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(50000000);
  gradient_structure::set_MAX_NVAR_OFFSET(1000);


DATA_SECTION
	
    int iseed;
    !!long int lseed=iseed;
    !!CLASS random_number_generator rng(iseed);
	
  // se leen los datos contenidos en *.dat. Se define si es matriz, vector,
   // entero o numero
  init_int nanos //numero de años 
  init_int ntallas  //nu,ero de tallas
  init_matrix indices(1,nanos,1,17)  //AQUI CONTENIDO DE vectore anuales
  //init_matrix Cl_cosecha(1,nanos,1,ntallas)  //Compocion de tallas de cosecha observada
  init_matrix Nlobs(1,nanos,1,ntallas)  //Composicion de tallas evaluacion
  //!!cout << "INDICES:"<< indices<<endl;exit(1);

  // ahora leo controles y opciones desde el archivo _opt.dat
  !! ad_comm::change_datafile_name("opmacha.ctl");
  matrix Wmed(1,nanos,1,ntallas)  //matriz de peso medio a la talla leyendo a y b
  vector arl(1,nanos)
  vector brl(1,nanos)
	  //!!cout << "brl"<<brl<<endl;exit(1); 
  init_vector Tallas(1,ntallas) //vector de tallas
  init_vector mat(1,ntallas)  //TODO: revisar talla de madurezmadurez a la talla
  init_vector Sel(1,ntallas) //Selectividad a la talla (talla > talla minima =1, 0 en otro caso)
	  
	  //!!cout<<"Sel "<< Sel << endl;exit(1);
  //Lee CV Desv R CV(K) CV(q)
  init_vector cvar(1,4)
  init_vector nmus(1,nanos) //para multinomial frecuencia de tallas
  init_vector dt1(1,nanos)  //Intervalos de tiempo dentro del año
  init_number dt2 //dt desove
  init_int    opt_qEval
  init_int    opt_Sel4 
  //Loo  k Lr sr Bp M h q
  init_vector parbiol(1,8)  //parametros biologicos
  init_int    opt_VB1  //opcion estima talla de reclutamiento (Lr)
  init_int    opt_VB2  //opcion estima desviacion talla de reclutamiento (sr)
  init_int    opt_VB3  //opcion estima K
  init_int    opt_VB4  //opcion estima Beta para la matriz de transicion
  init_int    opt_M    //opcion estima Mortalidad natural
	  //!!cout << "opt_M" << opt_M << endl;exit(1);
	//opcion Fase estimacion de parametros poblaciones
  init_int    opt_Rmed //fase estimacion de reclutamiento promedio   
  init_int    opt_devR //fase estimacion desvios de R
  init_int    opt_F    //fase estimacion de F
  init_int    nanos_proy //años futuros para proyeccion 
  init_int    npbr       //numero de pbrs
  init_vector mu(1,npbr) //Valor de tasa de explotacion para cada pbr
  init_number RPRmsy  //razon potencial reproductivo al MRS
  init_int    op_str  //Para utilizar rampa
  init_int SrType
  init_vector cambio(1,nanos_proy)

	  //!!cout<< "SrType " << SrType << endl; exit(1);

INITIALIZATION_SECTION
 // defino un valor inicial de log_reclutamiento promedio (factor de escala)
 log_Rmed     10.30
 // log_A50f_one 2.54
 // log_A50R_one 2.54
 // log_Df_one   0
 // dev_log_A50f 0
 // dev_log_Df   0
 //log_sf       -0.69
 // log_sr       0.37
 // log_Lr       2.07
 // log_bp       -2.41
  

PARAMETER_SECTION

//Se acota el dominio de cada parametro
// parametros selectividad Cosecha
 init_bounded_number log_A50R_one(1.6,4.1,opt_Sel4)  
 init_bounded_number log_DR_one(0,5,opt_Sel4)

// Recruits and F mortalities
// parametros reclutamientos y mortalidades)
 init_number log_Rmed(opt_Rmed)
 init_bounded_dev_vector log_desv_No(1,5,-5,5,opt_devR)
 init_bounded_dev_vector log_desv_Rt(1,nanos,-5,5,opt_devR)
 init_bounded_vector log_F(1,nanos,-6,1.39,opt_F) // log  mortalidad por pesca por flota
 init_bounded_number log_M(-2.3,0.69,opt_M)

// capturabilidades
 init_number log_qEval(opt_qEval)

// crecimiento
 init_bounded_number log_Lr(1.6,3,opt_VB1)
 init_bounded_number log_sr(-2.3,1,opt_VB2)
 init_bounded_number log_k(-1.6,0.18,opt_VB3)
 init_bounded_number log_bp(-2.3,-0.69,opt_VB4)

//Defino las variables de estado 
 vector yrs(1,nanos)
 vector Cobs(1,nanos)
 vector Cpred(1,nanos);
 vector Bobs(1,nanos)
 vector Bpred(1,nanos);
 
 vector cv1(1,nanos) //Coef variacion Biomasa
 vector cv2(1,nanos) //Coef variacion Capturas
 vector Unos_tallas(1,ntallas)
 vector Unos_yrs(1,nanos)
 //vector Reclutas(1,nanos)

 vector Bt(1,nanos);
 vector Bv(1,nanos);
 vector Bcru(1,nanos);
 vector prior(1,7);

 vector Neq(1,ntallas);
 vector likeval(1,10);
 vector SDo(1,nanos);
 vector Scru_1(1,ntallas);
 vector Scru_2(1,ntallas);
 vector log_A50f(1,nanos)

 vector log_Df(1,nanos)
 vector log_A50R2(1,nanos)

 vector log_DR2(1,nanos)
 
 vector Lesp(1,ntallas)
 vector sigmaL(1,ntallas)


 vector No(1,ntallas);
 matrix Sflo(1,nanos,1,ntallas)
 matrix Seval(1,nanos,1,ntallas)
 matrix F(1,nanos,1,ntallas)
 matrix Z(1,nanos,1,ntallas)
 matrix S(1,nanos,1,ntallas)
 matrix N(1,nanos,1,ntallas)
 matrix Nv(1,nanos,1,ntallas)
 vector Nop(1,ntallas) //  nueva lineas
 vector fyr(1,nanos)

 matrix Cl_pred(1,nanos,1,ntallas)
 matrix Nlpred(1,nanos,1,ntallas)
 matrix pobs(1,nanos,1,ntallas)
 matrix ppred(1,nanos,1,ntallas)
 matrix pobs_cru(1,nanos,1,ntallas)
 matrix ppred_cru(1,nanos,1,ntallas)
 matrix T(1,ntallas,1,ntallas)
 matrix Id(1,ntallas,1,ntallas);
 //nueva lineas
 number suma1
 number suma2
 number suma3
 number suma4
 number cuenta1
 number nm1

 number Linf
 number k
 number Lr
 number sr
 number M
 number qprior
 number Depl //


 vector pre(1,ntallas)
 vector delta(1,ntallas)
 vector alfa(1,ntallas)
 number bp
 vector dL(1,ntallas);
 vector Wm(1,ntallas); 

 sdreport_vector SD(1,nanos) // 
 //vector SA(1,nanos) //Biomasa adultos (desovantes al principio) 
 sdreport_vector RPR(1,nanos) // 
 sdreport_vector RPRlp(1,nanos)
 sdreport_vector Reclutas(1,nanos) // 
 //sdreport_vector Bsem(1,nanos)
 sdreport_number SSBo

 likeprof_number qacu 
 
  // los arreglos usados en la proyeccion
	 vector Np(1,ntallas)
	 vector Sp(1,ntallas)
	 vector Fp(1,ntallas)
	 vector Zp(1,ntallas)
		 
		 //MODELO OPRATIVO
		 vector Nlp(1,ntallas);
         vector Nlvp(1,ntallas);
         matrix Bp(1,npbr,1,nanos_proy) //Biomasa total
			 number Rp
		 matrix Bvp(1,npbr,1,nanos_proy)
		 matrix BDop(1,npbr,1,nanos_proy) //Biomasa desovante a comienzos de año
		 matrix BDp(1,npbr,1,nanos_proy) //Biomasa desovante despues de la mortalidad
 	     matrix Bcp(1,npbr,1,nanos_proy) //Biomasa total de la evaluación
		 matrix Bcvp(1,npbr,1,nanos_proy) //Biomasa explotable de la evaluacion
		 matrix CTP(1,npbr,1,nanos_proy) //Captura en numero desde la evaluacion directa
		 matrix YTP(1,npbr,1,nanos_proy) //Captura en peso desde la evaluación directa
		 matrix RPRp(1,npbr,1,nanos_proy)
		 matrix Fmortp(1,npbr,1,nanos_proy)
		 matrix Rproy(1,npbr,1,nanos_proy)

		 //matrix nproy(1,nanos_proy,1,ntallas)
			 //newton-raphson
		 number Kobs_captura
				 matrix temp(1,nanos,1,ntallas);
				 number R0
				 number S0
				 number alpha
				 number beta
				 number h
				 number SPR0
				 number ssb
				 vector Nspr(1,ntallas)
				 vector Rpred(1,nanos)
	 
  objective_function_value f


PRELIMINARY_CALCS_SECTION
// leo la matriz de indices

  yrs=column(indices,1);// asigno la 1 columna de indices años
  cv1=column(indices,5); //Coeficiente de variacion Bobs
  Bobs=column(indices,8); //Biomasa cruceros
  cv2=column(indices,16);  //Coeficiente de variacion Capturas
  Cobs=column(indices,17); //Capturas en peso
  arl=column(indices,12);
  brl=column(indices,13);

  Unos_yrs=1;// lo uso en operaciones matriciales con el año
  Unos_tallas=1;// lo uso en operaciones matriciales con el año
  //relacion longitud-peso
  for(int i=1;i<=nanos;i++)
  {

	  for(int l=1;l<=ntallas;l++)
	  {
		  Sflo(i,l) = Sel(l);
		  Wmed(i,l)=arl(i)*pow(Tallas(l),brl(i));
	  }
  }
  
  T=0;
  int j;
  
  for(j=1;j<=ntallas;j++){
	  Id(j,j)=1.0;
  }
  

RUNTIME_SECTION
  maximum_function_evaluations 200,1000,5000
  convergence_criteria  1e-3,1e-5,1e-6


PROCEDURE_SECTION


  Eval_Trans_talla_talla();
  Eval_selectividad();
  Eval_mortalidades();
  Eval_abundancia();
  Eval_capturas_predichas();
  Eval_deinteres();
  srr();
  Eval_logverosim();
  Eval_funcion_objetivo();
  //if(last_phase()){
	//  srr();
	  //Op_model();
	  //}
  if(mceval_phase()){
	  Op_model();}


FUNCTION srr

//Relacion stock-recluta

//Numero por recluta
 Nspr=pre*1*inv(Id-T*(exp(-1.*M)*Id));
 SPR0=sum(elem_prod(Nspr*exp(-dt2*M),elem_prod(Wmed(1),mat)));
 R0 = mean(Reclutas);//exp(log_Rmed);//-0.5*square(cvar(1)));
 S0 = R0*SPR0;
 //S0 = mean(SD);
 //R0 = S0/SPR0;
 if(SrType==1)
 {
   alpha = (S0/R0)*((1-h)/(4*h));
   beta = (5*h-1)/(4.*h*R0);
   for(int i=1;i<=nanos;i++)
   {
     Rpred(i)=SD(i)/(alpha + beta*SD(i));	
   }
 }
 else
 {
   alpha = (R0/S0)*exp(log(5*h)/0.8);
   beta = log(5*h)/(0.8*S0);
   for(int i=1;i<=nanos;i++)
   {
     Rpred(i)=alpha*SD(i)*exp(-beta*SD(i));	
   }
 }
 


FUNCTION Eval_Trans_talla_talla

  Linf=parbiol(1);
  k=parbiol(2);
  bp=parbiol(5);

  if(active(log_k)){k=mfexp(log_k);}
  if(active(log_bp)){bp=mfexp(log_bp);}

 int i, j;
 
// matriz de transicion modelo normal

  delta=(Linf-Tallas)*(1-mfexp(-k));// incremento en tallas
  Lesp=Tallas+delta; // talla esperada luego del crecimiento
  sigmaL=delta*bp;  

  for (i=1;i<=ntallas;i++){
     T(i)=exp(-0.5*square((Lesp(i)-Tallas)/sigmaL(i)));
     T(i)(1,i)=0;
  }

 
  for (j=1;j<=ntallas;j++){

    if(delta(j)<0){
    T(j)=0;
    T(j,j)=1;}
  } 

  T(ntallas,ntallas)=1;
  for (j=1;j<=ntallas;j++){
  T(j)/=sum(T(j));
  }
  h = parbiol(7);
  

FUNCTION Eval_selectividad
  int i;
  /*

  log_A50f=log_A50f_one+dev_log_A50f;
  log_Df=log_Df_one+dev_log_Df;
  cout<<log_A50f<<endl;
  cout<<log_Df<<endl;exit(1);
  // Selectividad logística flota
  for (i=1;i<=nanos;i++)
  {Sflo(i)=elem_div(Unos_tallas,(1+exp(-1.0*log(19)*(Tallas-exp(log_A50f(i)))/exp(log_Df(i)))));
   }
  // cout<<Sflo<<endl;exit(1);
  
  // parametros selectividad doble normal

  if (opt_Sel3>0){
   // selectividad doble_normal unica 
    for (i=1;i<=nanos;i++){
      Sflo(i)=mfexp(-1/(2*square(exp(log_sf(1))))*square((Tallas-exp(log_muf))));
      for (int j=1;j<=ntallas;j++){
       if(Tallas(j)>exp(log_muf)){
       Sflo(i,j)=mfexp(-1/(2*square(exp(log_sf(2))))*square(-1.*(Tallas(j)-exp(log_muf))));}
    }}}

	*/
// Selectividad cruceros

  Seval=1.0;
    
  if (opt_Sel4>0){
    for (i=1;i<=nanos;i++)
    {Seval(i)=elem_div(Unos_tallas,(1+exp(-1.0*log(19)*(Tallas-exp(log_A50R_one))/exp(log_DR_one))));
    }}
  
FUNCTION Eval_mortalidades
  int i;
  M=parbiol(6);
  if (active(log_M)){M=mfexp(log_M);}

  F=elem_prod(outer_prod(mfexp(log_F),Unos_tallas),Sflo);
  Z=F+M;
  S=mfexp(-1.0*Z);
  for(i=1;i<=nanos;i++){fyr(i)= max(F(i));}

FUNCTION Eval_abundancia
  int i, j;

  Lr=parbiol(3);
  sr=parbiol(4);

 if (active(log_Lr)){Lr=mfexp(log_Lr);}
 if (active(log_sr)){sr=mfexp(log_sr);}


// genero la composicion de tallas del reclutamiento
  pre=exp(-0.5*square((Tallas-Lr)/sr));
  pre/=sum(pre);

// genero una estructura inicial de equilibrio en torno a Z del primer año;
  Reclutas=mfexp(log_Rmed+log_desv_Rt-0.5*square(cvar(1)));
  Neq=pre*Reclutas(1);

  for (j=1;j<=5;j++)
  {Neq=elem_prod(Neq,exp(-1.*Z(1)))*T+pre*exp(log_Rmed+log_desv_No(j)-0.5*square(cvar(1)));}

// genero la poblacion en equilibrio virginal de LP;
  //No(1)=pre*mfexp(log_Rmed)*mfexp(-0.5*square(cvar(1)));
  //for (j=2;j<=10;j++)
  //{No(j)=(No(j-1)*exp(-1.*M))*T;}
    //No=pre*mfexp(log_Rmed-0.5*square(cvar(1)))*inv(Id-T*(exp(-1.*M)*Id));
    No=pre*Reclutas(1)*inv(Id-T*(exp(-1.*M)*Id));
// dinamica anual

  N(1)=Neq;// primer año

  for (i=2;i<=nanos;i++)
  {N(i)=elem_prod(N(i-1),S(i-1))*T+pre*Reclutas(i);} //
//  {N(i)=elem_prod(N(i-1)*T,S(i-1))+pre*Reclutas(i);} //


FUNCTION Eval_capturas_predichas
	for(int i=1;i<=nanos;i++)
      {
			  temp(i)=dt1(i);
       }
 //temp=Unos_yrs*dt1;
// matrices de capturas predichas por talla y año
  Cl_pred=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));

// matrices de cruceros predichas por talla y año
  Nlpred=elem_prod(elem_prod(N,exp(-elem_prod(temp,Z))),Seval); //TODO:dt1 es un vector Verificar temp y outer_prod!	
  //cout << "Nlpred " << Nlpred << endl;exit(1);
// Composicion por tallas Evaluacion directa y estimada
  pobs_cru=elem_div(Nlobs,outer_prod(rowsum(Nlobs),Unos_tallas)+1e-5);
  ppred_cru=elem_div(Nlpred,outer_prod(rowsum(Nlpred),Unos_tallas)+1e-5);

// vectores de desembarques predichos por año
  Cpred=rowsum(elem_prod(Cl_pred,Wmed));



FUNCTION Eval_deinteres


// Rutina para calcular RPR
  Nv(1)=N(1);// solo para empezar los calculos
  
 for (int i=2;i<=nanos;i++)
  {Nv(i)=(Nv(i-1)*mfexp(-1*M))*T+pre*Reclutas(i);}

  Bt=rowsum(elem_prod(N,Wmed));
  Bv=rowsum(elem_prod(elem_prod(N,Sflo),Wmed));
  Bpred=rowsum(elem_prod(Nlpred,Wmed));
  
  qacu=parbiol(8);

  if (active(log_qEval)){qacu=mfexp(log_qEval);}
  
  Bpred=qacu*Bpred;

  //BIOMASA DESOVANTE
  SD=rowsum(elem_prod(elem_prod(N,exp(-dt2*Z)),elem_prod(Wmed,outer_prod(Unos_yrs,mat))));
  SDo=rowsum(elem_prod(Nv*exp(-dt2*M),elem_prod(Wmed,outer_prod(Unos_yrs,mat))));
  SSBo=sum(elem_prod(No*exp(-dt2*M),elem_prod(Wmed(1),mat)));
  //SA=rowsum(elem_prod(N,outer_prod(Unos_yrs,elem_prod(Wmed,outer_prod(Unos_yrs,mat)))));

 /*
  SD=rowsum(elem_prod(N,outer_prod(Unos_sem,elem_prod(Wmed,msex))));
  SDo=rowsum(elem_prod(Nv,outer_prod(Unos_sem,elem_prod(Wmed,msex))));
  SSBo=sum(elem_prod(colsum(No),elem_prod(Wmed,msex)));
 */


  RPR=elem_div(SD,SDo);
  RPRlp=SD/SSBo;


FUNCTION Eval_logverosim
// esta funcion evalua el nucleo de las -log-verosimilitudes marginales para
// series con datos 0.
  int i;

  suma1=0; suma2=0; 

  for (i=1;i<=nanos;i++)
  {
   if (Bobs(i)>0){
    suma1+=square((log(Bobs(i))-log(Bpred(i)))/cv1(i));}
	if (Cobs(i)>0){
	suma2+=square((log(Cobs(i))-log(Cpred(i)))/cv2(i));}
  }


FUNCTION Eval_funcion_objetivo

  int i,j;

 // se calcula la F.O. como la suma de las -logver
 // lognormalgraf
  likeval(1)=0.5*suma1;//Cruceros
  //likeval(2)=0.5*suma2;//CPUE
  likeval(3)=0.5*suma2;//CTP  

 // multinomial flota
 likeval(4)=0;
 //suma1=0;
 for(i=1;i<=nanos;i++)
 {
	 suma1=sum(elem_prod(pobs_cru(i),log(ppred_cru(i))));
     likeval(4)+=-1.*nmus(i)*suma1; 	
 } 

// multinomial cruceros
  //likeval(5)=-1.*nmus(2)*sum(elem_prod(pobs_cru,log(ppred_cru)));

// Priors
  prior(1)=0.5*norm2(log_desv_Rt)/square(cvar(1));
  prior(3)=0.5*norm2(log(Reclutas)-log(Rpred))/square(cvar(1));
  if (active(log_k)){// si estima k
  prior(2)=0.5*square((log_k-log(parbiol(2)))/cvar(2));}

  if (active(log_qEval)){
  prior(4)=0.5*square((exp(log_qEval)-parbiol(8))/cvar(3));}

  //if (active(dev_log_A50f)){
  //prior(6)=0.5*norm2(dev_log_A50f)/square(cvar(4));}

   f=sum(likeval)+sum(prior);

FUNCTION Op_model
	
	dvector dev_rec_fut(1,nanos_proy);
    dvector dev_surv_fut(1,nanos_proy);

	dev_rec_fut.fill_randn(rng);
	dev_surv_fut.fill_randn(rng);
	dvariable sigma_cru=0.1;	
	
   for(int j=1;j<=npbr;j++)
   {
	     Np=N(nanos);
		 Sp=S(nanos);
		 Zp=Z(nanos);
		 Nop=Nv(nanos);
		 ssb=SD(nanos);
		 Wm = Wmed(1);
		 //cout << "Np " << Np << endl;
		 //cout << "Sp" << Sp << endl;
		 for(int i=1;i<=nanos_proy;i++)
		 {
			 if(SrType==1)
			 {
				 Rp=(ssb/(alpha+beta*ssb))*mfexp(dev_rec_fut(i)*cvar(1)*cambio(i)-0.5*square(cvar(1))); //mfexp(dev_rec_fut(i)*cvar(1)*cambio(i)-0.5*square(cvar(1)));			 	
			 }
			 else
			 {		  Rp=alpha*ssb*exp(-beta*ssb)*mfexp(dev_rec_fut(i)*cvar(1)*cambio(i)-0.5*square(cvar(1))); //mfexp(dev_rec_fut(i)*cvar(1)*cambio(i)-0.5*square(cvar(1)));			 	
			 }
			 //cout<< "Rp " << Rp << endl;
			 Rproy(j,i)=Rp;
			 Np=elem_prod(Np,Sp)*T + pre*Rp;  
			 //cout << "Np" << Np << endl;
			 Nop=Nop*mfexp(-1.0*M)*T + pre*Rp;
			 //abundancia poblacional por talla disponible evaluación a comienzos de año
			 Nlp =elem_prod(Np,Seval(nanos));
			 //abundancia por talla explotable o vulnerable
			 Nlvp= elem_prod(Nlp,Sflo(nanos)); //Se usa Nlp 
			 Bp(j,i)=sum(elem_prod(Np,Wm)); //biomasa total de machas
			 //cout << "Wmed" << Wmed(nanos) << endl;
			 //cout << "Bp" << Bp << endl;
			 
			 Bvp(j,i) = sum(elem_prod(Nlvp,Wm)); //biomasa explotableo vulnerable real
			 Bcp(j,i)=sum(elem_prod(Nlp,Wm)); //biomasa total evaluada
			 Bcp(j,i)=qacu*Bcp(j,i)*mfexp(dev_surv_fut(i)*sigma_cru);
			 Bcvp(j,i)=sum(elem_prod(Nlvp,Wm)); //Biomasa explotable para calculo de CTP
			 Bcvp(j,i)=qacu*Bcvp(j,i)*mfexp(dev_surv_fut(i)*sigma_cru);
			 CTP(j,i)=mu(j)*qacu*sum(Nlvp)*mfexp(dev_surv_fut(i)*sigma_cru); //Captura  en numero
			 YTP(j,i)=mu(j)*Bcvp(j,i); //Captura en peso empirica mu es la tasa de explotacion
			 //Newton-Raphson
			 Kobs_captura=YTP(j,i);
			 dvariable ffpen=0.0;
			 dvariable SK=posfun((Bvp(j,i)-YTP(j,i))/Bvp(j,i),0.01,ffpen);
			 Kobs_captura=Bvp(j,i)-SK*Bvp(j,i);
			 do_Newton_Raphson_for_mortality(j,i);
			 //Fp=Fnew*Sflo(nanos);
			 Zp=M+Fp;
			 Sp=mfexp(-1*Zp);
			 //cout << "Sp " << Sp << endl;exit(1);
			 Fmortp(j,i)=max(Fp);
			 BDp(j,i) = sum(elem_prod(elem_prod(Np,exp(-dt2*Zp)),elem_prod(mat,Wm)));
 			 BDop(j,i)= sum(elem_prod(Nop*mfexp(-dt2*M),elem_prod(mat,Wm)));   //  nueva lineas
			 ssb=BDp(j,i); 
 			 RPRp(j,i)=BDp(j,i)/BDop(j,i);
			 //cout << "BDp " << BDp << endl;exit(1);
		 }
   }
   
   rep1 << Bp << endl; //Biomasa poblacional de macha real
   rep2 << Bvp << endl; //Biomasa poblacional disponible para la evaluacion real
   rep3 << Bcp << endl; //Biomasa evaluada por metodo directo estimador de Bp
   rep4 << Bcvp << endl; //Biomasa explotable fraccion de Bcp para calculo de CTP
   rep5 << BDp << endl; //Biomasa desovante real de macha
   rep6 << RPRp << endl; //razón de potencial reproductivo anual
   rep7 << CTP << endl; //captura en numero
   rep8 << YTP << endl; //captura en peso
   rep9 << Fmortp << endl;   //mortal
   rep10 << SD << endl;
   rep11 << fyr << endl;
   rep12 << Reclutas << endl;
   rep13 << Rproy << endl;
   

FUNCTION void do_Newton_Raphson_for_mortality(int j, int i)

     dvariable Fold = Kobs_captura/Bvp(j,i);
     dvariable Fnew ;
     for (int ii=1;ii<=5;ii++)
     {
         dvariable ZZ = Fold + M;
         dvariable XX = exp(-1.*ZZ);
         dvariable AA = Fold * (1. - XX);
         dvariable BB = ZZ;
         dvariable CC = 1. + (Fold - 1) * XX;
         dvariable dd = 1.;
         dvariable FX = AA / BB - Kobs_captura/Bvp(j,i);
         dvariable FPX = (BB * CC - AA * dd) / (BB * BB);
         Fnew = Fold - FX / FPX;
         Fold = Fnew;
     }
     Fp=Fnew*Sflo(nanos);

   
   /*
  FUNCTION Eval_CTP


  // CTP año en curso

 for (int j=1;j<=npbr;j++){

  Fp=Sflo(nanos)*Fpbr(j);//
  Zp=Fp+M;
  Sp=mfexp(-1.*Zp);
  CTPac=elem_prod(elem_div(Fp,Zp),elem_prod(Np,(1-exp(-Zp))));
  YTPac(j)=sum(elem_prod(CTPac,Wmed));
  }
  

 for(int j=1;j<=npbr;j++)
 {
	 Np=N(nanos);
	 Sp=S(nanos);
	 Zp=Z(nanos);
	 Nv2=Nv(nanos); //  nueva lineas
	 Depl=RPR(nanos);
	 for(int i=1;i<=nanos_proy;i++)
		{
			Np=elem_prod(Np,Sp)*T+pre*mfexp(log_Rmed);
			Nv2=Nv2*mfexp(-1.0*M)*T+pre*mfexp(log_Rmed);//  nueva lineas
			Fp=Sflo(nanos)*Fpbr(j);//
			Fproy(i,j)=Fpbr(j);
			if(op_str==0)
			{
				if(RPRp(i,j)<RPRmsy)
				{
					if(Depl<RPRmsy)
					{
						Fproy(i,j)=Fpbr(j)*Depl/RPRmsy;
						Fp=Sflo(nanos)*Fproy(i,j);//
						respue(i,j)=1;
					}
				}
			Zp=Fp+M;
			Sp=mfexp(-1.*Zp);
			BPA(i,j)=sum(elem_prod(Np,(elem_prod(Wmed,outer_prod(Unos_yrs,mat)))));
			Bp(i,j)=sum(elem_prod(elem_prod(Np,exp(-dt2*Zp)),elem_prod(Wmed,outer_prod(Unos_yrs,mat))));
			Bvp(i,j)= sum(elem_prod(elem_prod(Nv2*mfexp(-dt2*M),mat),Wmed));   //  nueva lineas
			RPRp(i,j)=Bp(i,j)/Bvp(i,j);
			Depl=RPRp(i,j);
			CTP=elem_prod(elem_div(Fp,Zp),elem_prod(Np,(1-Sp)));
			YTP(i,j)=sum(elem_prod(CTP,Wmed));
			nproy(i)=Np;
		}
	}
  ofstream tac_out("ctps.dat");
  tac_out<< 
  YTP <<
  endl;
  tac_out.close();
 
  ofstream grad_out("grad_final.dat"); 
  grad_out<< objective_function_value::pobjfun->gmax <<endl;
  grad_out.close();
  */

REPORT_SECTION

// si quiero que este proceso lo haga despues de la hessiana, no tendre IC
//  Bsem=rowsum((elem_prod(N*Prob_talla,outer_prod(Unos_sem,Wmed))));
   suma1=0; suma2=0;nm1=1;cuenta1=0;

  for(int i=1;i<=nanos;i++)
  {
	  if(sum(pobs_cru(i))>0)
	  {
         suma1=sum(elem_prod(ppred_cru(i),1-ppred_cru(i)));
		 suma2=norm2(pobs_cru(i)-ppred_cru(i));
		 nm1=nm1*suma1/suma2;
		 cuenta1+=1;
	}}
	nm1=pow(nm1,1/cuenta1);

  rep(nm1);
  rep(Bobs);
  rep(Bpred);
  rep(Cobs);
  rep(Cpred);
  rep(Reclutas);
  rep(log_desv_Rt);
  rep(Bt);
  rep(Bv);
  rep(SD);
  rep(RPR);
  rep(RPRlp);
  rep(Tallas);
  rep(pre);
  rep(S);
  rep(pobs_cru);
  rep(ppred_cru);
  rep(N);
  rep(F);
  rep(Sflo);
  rep(Seval);
  rep(Linf);
  rep(k);
  rep(Lr);
  rep(sr);
  rep(bp);
  rep(qacu);
  rep(M);
  rep(exp(log_Rmed));
  rep(T);
  rep(Wmed);
  rep(No);
  rep(SSBo);
  rep(temp);
  rep(log_desv_Rt)
  rep(S0);
  rep(R0);
  rep(SPR0)
  rep(h)
  rep(alpha)
  rep(beta)
  rep(Nspr)
  rep(likeval)
  rep(f)
  rep(prior);
  rep(mat);
  rep(fyr);
  rep(cambio);
  

GLOBALS_SECTION
    #include <admodel.h>
    #include <fenv.h>
    #undef rep
    #define rep(object) report << #object "\n" << object << endl;
	
	ofstream rep1("01Btota_mcmc.txt",ios::app);
    ofstream rep2("02Bvuln_mcmc.txt",ios::app);
	ofstream rep3("03Bcru_mcmc.txt",ios::app);
	ofstream rep4("04BcruVuln_mcmc.txt",ios::app);
	ofstream rep5("05Bdes_mcmc.txt",ios::app);
	ofstream rep6("06RPR_mcmc.txt",ios::app);
	ofstream rep7("07CTP_mcmc.txt",ios::app);
	ofstream rep8("08YTP_mcmc.txt",ios::app);
	ofstream rep9("09Fmort_mcmc.txt",ios::app);
	ofstream rep10("10SD_historia.txt",ios::app);
	ofstream rep11("11Fmort_historia.txt",ios::app);
	ofstream rep12("12Reclutas_historia.txt",ios::app);	
    ofstream rep13("13Recluta_mcmc.txt",ios::app);

