#### CONTROL DE DATOS BIOLOGICOS Y DE AJUSTE
#grupos	de	Tallas
5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90	
#	pesos	medios	a	la	talla estimados internamente
#	ojiva	de	madurez	sexual
0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.5 0.9 1.0 1.0 1.0 1.0 1.0 1.0 1.0
# Selectividad de la pesqueria
0  0  0  0  0  0  0  0  0  1  1  1  1  1 1 1 1 1
# Coeficientes de variacion y tamaños de muestra
#---------------------------------------------------
#    cv_devRmed       cv_k       cv_q_cru    cv_desv_L50
     0.5           0.1         0.01          0.15
#  (nmus) Tamaños de muestra composiciones de tallas evaluacion directa
400 50 180 180 0 400 400s
#---------------------------------------------------
#  (dt) dt del año de las evaluaciones y del desove
# 2011 2012  2013    2014 2015 2016   2017
0.5	   0.25	 0.333	0.083	 0    0.083   0.083
# dt desove
0.81
# Fases estimacion capturabilidad
#---------------------------------------------------
#   (opt_q <0 es qprior) capturabilidad del Evaluacion directa
4
# Fases/Opciones de selectividad
#---------------------------------------------------
#  (opt_Sel) Selectividad de cruceros (si <0, Scru=1)
-2
#  Parámetros biológicos
#---------------------------------------------------
#  Loo    k      Lr    sr     Bp      M	         h     q   q_eval(prior)																					
#  93.4   0.20   20.0    1.5   0.2    0.3   0.7    0.9    #Rubilar et al. (2001)
  93.4   0.25   20.0    1.5   0.2    0.3   0.7    0.99    #Rubilar et al. (2001)
#  Fase de estimacion de parámetros crecimiento
#---------------------------------------------------
# (opt_VB1) Estimacion de Lr (es la primera talla modal)																						
3
# (opt_VB2) Estimacion de desviacion recluta (sr)
3
# (opt_VB3) Estimacion de k
3
# (opt_VB4) Estimacion de Bp
3
# (opt_M) Estimación de  M																						
-3																						
#  Fase de estimacion de parámetros poblacionales
#---------------------------------------------------
# (opt_Rmed) Estimacion de Rmed																						
1																						
# (opt_devR) Estimacion desvios R																				
2																						
# (opt_F) Estimación de F																						
2																						
#  Opciones de proyeccion
#---------------------------------------------------
# N años futuro
20   
# n de pbrs
6
# Tasa de explotacion
0 0.1  0.15 0.2  0.25 0.3 
#RPRmsy
#0.60
0.45
#para usar rampa =0
1
#TypeSR (1: Bev-Holt; 2: Ricker)
1
#switch
1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1


