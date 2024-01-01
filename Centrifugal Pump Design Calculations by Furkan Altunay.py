#Santrifuj Pompa Hesap by Furkan Altunay
import math
from math import sqrt
from math import pow
from sympy import symbols, solve, Eq


# n (d/dk)
n=1450
# Q52(lt/sn)
Q52=65
Q=Q52*0.001
# manometrik yukseklik Hm (m)
Hm=23
g=9.81
hidroeffic=0.85

# a) ozgul hiz
ns=round(3.65*n*((Q**(1/2))/(Hm**(3/4))), 2)
print ("ozgul Hiz ns (d/dak): ", ns)

# b) Pompa mil gucu ve cevirici motor gucu
effic52=(64.9 + 0.157*ns + -4.37E-04*(ns**2) + 2.67E-07*(ns**3))*0.01
effic=effic52-(effic52*(5/100))
print("effic (%): ",effic)
Ne=round((((10**3)*Q*Hm))/(75*effic), 1)
print ("Mil gucu (efektif guc) (Beygir gucu): ", Ne)
#Emniyet katsayisi 1.15 alinarak, elektrik motorunun gucu
EK=1.15
Nm=EK*Ne
print("Elektrik motorunun gucu (B.G.): ", Nm)

# c) Mil capi: kayma gerilmesi karbonlu celik icin 120 kg/cm2 secilerek
dmil = round(14.4 * math.pow(Ne / n, 1/3))
dmil_check = 14.4 * math.pow(Ne / n, 1/3)
print(dmil_check)
print("Mil capi dmil (cm): ", dmil)
dmil2=dmil*10 # mm 

#Gobek capi
dgobek=1.5*dmil2
print("Gobek capi dgobek (mm): ", dgobek)

# d) Pompa carki giris boyutlari
#Emme borusundaki su hizi
Ve=1.07 + 0.0416*Q52 + -1.83E-04*(Q52**2) + 3.87E-07*(Q52**3)
print("Ve (m/s): ", Ve)

#Emme borusu capi De
De=math.sqrt((4*Q)/(math.pi*Ve))
print("Emme Borusu capi De (m): ", De)
#Bu deger en yakin norm cap degerine getirilerek yuvarlanma yapilir
De2=round(De, 1)*1000
print("Emme Borusu capi Yuvarlanmis De (mm): ", De2)

Ve2=round((4*Q)/(math.pi*(De2/1000)**2), 2)
print("Borudaki Ve hizi (m/s): ",Ve2)

#D0 cark giris agzi capi, De capina esit alinabilir
D0=De2
print("cark Giris Agzi capi De (mm) :", D0)

#D0 capindaki kesitte C0 hizi, cark icindeki kacaklar goz onune alinarak bulunacak QT toplam debi yardimiyla ve sureklilik denkleminden hesaplanabilir. 

effic2=0.96 #Kacak verim secildi

QT=round((Q/effic2)*1000)
print("Toplam Debi QT (lt/s): ", QT)

C0=round((4*(QT/1000))/(math.pi*((D0/1000)**2-(dgobek/1000)**2)), 2)
print("D0 capindaki kesitte C0 hizi (m/s): ", C0)

#Kanatlardan dolayi kesit kuculeceginden D1 cark giris ortalama capina ait kesitteki C1 mutlak hizi, C0 hizina gore bir miktar arttirilarak asagidaki deger bulunur.
C1=round(1.1*C0, 2)
print("C1 mutlak hizi (m/s): ", C1)

Cm1=C1 #Dik giris icin Cm1=C1 olur. 
print("Giris meridyen hizi (m/s), Cm1: ", Cm1)

kem1=0.0659 + 1.15E-03*ns + -2.77E-06*(ns**2) + 3.91E-09*(ns**2)
Cm1_check=round(kem1*(math.sqrt(2*9.81*Hm)), 2)
print("Cm1_check (m/s): ", Cm1_check)

#D1 cark giris ortalama hizi o- Schulz katsayisi 0.92 secilerek
o=0.92

D1=o*De2/1000
print("cark giris ortalama capi (m), D1: ", D1)
D1d=D0+4
print("Kanat Giris Dis capi (mm), D1d: ", D1d)
D1i=2*D1*1000-D1d
print("Kanat Giris Ic capi (mm), D1i: ", D1i)

#beta1 kanat girisinde, akiskana ait giris acisi (optimum calisma noktasinda on donmesiz giris oldugu kabul edilerek alfa=90 derece, Cu1=0 ve C1=Cm1) D1 giris ortalama capina gore cevresel hiz U1 ve Cm1 giris meridyen hizi yardimiyla bulunur. 
U1= round((math.pi*D1*n)/(60), 2)
print("cevresel hiz (m/s), U1: ", U1)

beta1_old=math.atan(Cm1/U1)
ondalik_kisim = int(beta1_old* 1000) % 1000  
beta11=(ondalik_kisim)
beta1=beta11*0.001
print("Kanat Girisinde Akiskana Ait Giris Acisi, beta1: ", beta1)

aci_derece = beta1 * (180 / math.pi)  
derece_tam = int(aci_derece)  
kalan_derece = aci_derece - derece_tam  
kalan_dakika = kalan_derece * 60 
kalan_dakika_tam = round(kalan_dakika)  
print("beta1: "f"{derece_tam} derece {kalan_dakika_tam} dakika")

#Kaskad kanat etkisi dolayisiyla akiskanin bu aci ile carka girebilmesi icin kanat giris acisi delta kadar buyutulmelidir. delta 2 derece 23 dakika olarak secilerek
delta_derece = 2
delta_dakika = 23
delta_radyan = (delta_derece + delta_dakika / 60) * (math.pi / 180)
beta1k=beta1+delta_radyan
aci_derece2 = beta1k * (180 / math.pi)  
derece_tam2 = int(aci_derece2) 
kalan_derece2 = aci_derece2 - derece_tam2 

if kalan_derece2 >= 0.5:
    derece_tam2 += 1
    kalan_dakika_tam2 = 0
else:
    kalan_dakika_tam2 = round(kalan_derece2 * 60)

if kalan_dakika_tam2 == 60:
    derece_tam2 += 1
    kalan_dakika_tam2 = 0

print(f"beta1k: {derece_tam2} derece {kalan_dakika_tam2} dakika")

#Kanat giris kenari bakimindan beta1k acisi kontrol edilirse,

U1i=round((math.pi*(D1i/1000)*n)/60, 2)
U1d=round((math.pi*(D1d/1000)*n)/60, 2)
print("U1i (m/s): ", U1i)
print("U1d (m/s): ", U1d)

beta1i=round(Cm1/U1i,3)
beta1d=round(Cm1/U1d, 3)
print("beta1i (radyan): ", beta1i)
print("beta1d (radyan): ", beta1d)
# Radyan cinsinden aciya cevirme
beta1i_aci = math.degrees(math.atan(beta1i))
beta1d_aci = math.degrees(math.atan(beta1d))

# Acilari derece ve dakika olarak ayirma
beta1i_derece = int(beta1i_aci)
beta1i_dakika = round((beta1i_aci - beta1i_derece) * 60)

beta1d_derece = int(beta1d_aci)
beta1d_dakika = round((beta1d_aci - beta1d_derece) * 60)

print("beta1i: ", beta1i_derece, "derece", beta1i_dakika, "dakika")
print("beta1d: ", beta1d_derece, "derece", beta1d_dakika, "dakika")

# Aci kontrolu
if abs(beta1i - beta1d) < 6:
    print("Kanat silindirik olarak imal edilebilir.")
else:
    print("Kanat silindirik olarak imal edilemez.")

# b1 cark ginis genisligi
#Kanat kalinligi e ve kanat sayisi Z degerleri henuz belirlenmemis oldugu icin lambda1 degeri ilk yaklasimda 0.7 olarak secilirse, sureklilik denklemi yardimiyla b1 degeri yaklasik olarak bulunabilir.
lambda1=0.7
b1=round((QT/1000)/(math.pi*D1*lambda1*Cm1), 2)
print("cark Giris Genisligi Tahmini (m), b1_tahmin: ", b1)  
#Bu deger, cark cikisina ait boyutlar hesaplandiktan sonra kontrol edilecektir. 

# e) carkin cikis kosullari ve boyutlarinin saptanmasi
#Basinc katsayisi
pressurecoeff=1.47 + -6.25E-03*ns + 1.1E-05*(ns**2) + -2.57E-09*(ns**3)

U2=round(math.sqrt((2*g*Hm)/pressurecoeff), 2)
print("cikis cevresel Hizi (m/s), U2: ",U2)

D2=round((60*U2)/(math.pi*n), 3)
print("cark cikis capi (m), D2: ", D2)

ku2=0.973 + -6.12E-04*ns + 1.07E-05*(ns**2) + -1.72E-08*(ns**3)
U2_new=round(ku2*(math.sqrt(2*g*Hm)), 2)
print("cikis cevresel Hizi (m/s), U2_new: ",U2_new)

D2_new=round((60*U2_new)/(math.pi*n), 3)
print("cark cikis capi (m), D2_new: ", D2_new)

D2_new2_text = input("cark cikis capini Kac Almak Istiyorsunuz? (m) :")
D2_new2 = float(D2_new2_text)

U2_new2=round((D2_new2*math.pi*n)/60, 2)
print("cikis cevresel Hizi (m/s), U2_new2: ",U2_new2)

Drate=round(D1/D2_new2, 3)
print("D1/D2 oraninin sinir degeri: ", Drate)

kmc2=0.0527 + 7.34E-04*ns + -8.15E-07*(ns**2) + 7.22E-10*(ns**3)
Cm2=round(0.124*math.sqrt(2*g*Hm), 2)
print("cikis Meridyen Hizi (m/s), Cm2: ", Cm2)


Cu2=round(((g*Hm)/(U2_new2*hidroeffic)), 2)
print("Tegetsel Hiz Bilesenleri (m/s), Cu2: ", Cu2)

beta2_radyan=round(math.atan(Cm2/(U2_new2-Cu2)), 4)
beta2_derece = round(math.degrees(beta2_radyan),0)
print("beta2 (radyan):", beta2_radyan)
print("Akiskana Ait carktan cikis Acisi (derece), beta2:", beta2_derece)


# beta2k kanat cikis acisi: (Kaskat kanat etkisiyle akiskanin beta2 acisi ile carki terkedebilmesini saglamak uzere kanada verilmesi gerekli olan aci)
beta2k=20
print("Kanat cikis Acisi Kabulu, beta2k: ", beta2k)
beta2k_radyan = round(math.radians(beta2k), 3)
beta1k_v2= round(aci_derece2, 2)

tolerance = 0.5  # Istenen hassasiyet seviyesini ayarlayabilirsiniz
max_iterations = 100000
for iteration in range(max_iterations):
    # Beta2k degerini radyan cinsinden hesaplayin
    beta2k_radyan = math.radians(beta2k)
    
    Z = round(6.5 * ((D2_new2*1000+D1*1000)/(D2_new2*1000-D1*1000)) * math.sin(math.radians((beta1k_v2+beta2k)/2)), 2)
    Z_new=round(Z, 0)
    K=round(1+(((1.2*(1+beta2k_radyan))/(Z_new))*(1/(1-(D1/D2_new2)**2))), 3)
    Cu2sonsuz = round((K * Cu2), 3)

    # Yeni beta2k degerini hesaplayin
    beta2k_radyan2 = math.atan(Cm2 / (U2_new2 - Cu2sonsuz))
    beta2k_derece2 = math.degrees(beta2k_radyan2)

    # Aci kontrolu yapin ve istenen hassasiyete ulasildiysa donguyu sonlandirin
    if abs(beta2k - beta2k_derece2) < tolerance:
        print("beta2k Acisini Secmek Uygundur ---->", round(beta2k, 2))
        break
    else:
        beta2k = beta2k_derece2  # Yeni beta2k degeriyle donguyu tekrarlayin

# Maksimum iterasyon sayisina ulasildiginda veya istenen hassasiyete ulasilamadiginda uygun bir mesaj yazdirin
if iteration == max_iterations - 1:
    print("Hassasiyet seviyesine ulasilamadi. Daha fazla iterasyon veya baska bir yaklasim gerekebilir.")

print("Kanat Sayisi, Z: ", Z_new)
print("K degeri: ", K)
#Cu2sonsuz, sonsuz kanat halindeki tegetsel hiz bileseni
print("Sonsuz Kanat Halindeki Tegetsel Hiz (m/s), Cu2sonsuz: ", Cu2sonsuz)

# b2 cikis genisliginin hesabi
# Kanat kalinligoi e=5 mm alinirsa, lambda2 kanat daralma katsayisi
e=5
print("Kanat Kalinligi Alindi (mm) --->", e)
lambda2= round(1-((Z*e/(beta2k_radyan))/(math.pi*D2_new2*1000)), 2)
print("lambda2: ",lambda2)

b2=round((QT/1000)/(math.pi*D2_new2*lambda2*Cm2), 4)
print("cikis Genisligi (m), b2: ", b2)

b2_imalat=round(b2*1000+1, 0 )
print("cikis Genisligi Imalat Kolayligi acisindan alinmali(mm), b2: ", b2_imalat)

lambda1_gercek= round(1-((Z*e/(beta1k))/(math.pi*D1*1000)), 3)
print("lambda1_gercek: ",lambda1_gercek)

b1_gercek=round((QT/1000)/(math.pi*D1*lambda1_gercek*Cm1), 3)
print("Giris Genisligi (m), b1_gercek: ", b1_gercek)

#(Z.L)min kontrolu
betasonsuz=round((2*Cm2/(U2_new2+U1-Cu2)), 3)
betasonsuz_radyan=(math.atan(betasonsuz))
betasonsuz_derece = round(math.degrees(betasonsuz_radyan),2)
derece = int(betasonsuz_derece)
dakika = round((betasonsuz_derece - derece) * 60, 0)
print("betasonsuz (radyan): ", betasonsuz_radyan)
print("betasonsuz: ", f"{derece} derece {dakika} dakika")

Wsonsuz=(U1+U2_new2-Cu2)/(2*math.cos(betasonsuz_radyan))
print("Wsonsuz (m/s): ",Wsonsuz)

betam=(beta1k_v2+beta2k)/2
print("betam (derece): ", betam)
betam_radyan = math.radians(betam)

l=round((D2_new2-D1)/(2*math.sin(betam_radyan)), 3)
print("l (m): ", l)

t_l=round((math.pi*(D2_new2+D1))/(2*Z_new*l), 3)
print("t/l: " ,t_l)

#t/l ve betam degerleri icin S'p/Sp degeri bulunur.
S=1.51 
Z_l_min=round((2*math.pi*D2_new2*Cu2)/(1.5*Wsonsuz*S), 3)
print("Z.l min: ",Z_l_min)

Z_l= round(Z_new*l, 3)
print("Z.l: ",Z_l)
# Z.L min kontrolu
if Z_l_min < Z_l:
    print("(Z.l)min kontrolunden gecti")
else:
    print("(Z.l)min kontrolunden gecemedi")

#Difuzorun gerekli olup olmadiginin kontrolu
alfa2=round((Cm2/Cu2), 3)
alfa2_radyan=(math.atan(alfa2))
alfa2_derece = round(math.degrees(alfa2_radyan),2)
derece2 = int(alfa2_derece)
dakika2 = round((alfa2_derece - derece2) * 60, 3)
print("alfa2: ", f"{derece2} derece {dakika2} dakika")

if 10 < derece2:
    print("Difuzore Gerek Yok!")
else:
    print("Difuzore Gerek Var!")

#Salyangoz Hesabi
D3=D2_new2*1000+9
print("D3 (mm): ",D3) 
t=0.013
sb=8 #salyangoz bolunme sayisi

A=Cu2*(D2_new2/2)
print("Alanlar Sabitesi, A (m2/s)",A)

Qx1=Q52/sb/1000
R3=(D3/2)/1000
print("R3 (m): ",R3)

rx1 = symbols('rx1')
equation = Eq(A * math.pi * rx1**2 - Qx1 * rx1 - Qx1 * (R3 + t), 0)
rx1_solutions = solve(equation, rx1)
# Pozitif cozumu bulun
positive_rx1 = [sol.evalf() for sol in rx1_solutions if sol > 0]
rx1=round(positive_rx1[0], 3)
print("rx1 (m): ", rx1)

d1=rx1*2*1000
print("d1 (mm): ", d1)

print("1. kesit Qx1 (m3/s): ", Qx1)
Qx2=2*Q52/8/1000
print("2. kesit Qx2 (m3/s): ", Qx2)
Qx3=3*Q52/8/1000
print("3. kesit Qx3 (m3/s): ", Qx3)
Qx4=4*Q52/8/1000
print("4. kesit Qx4 (m3/s): ", Qx4)
Qx5=5*Q52/8/1000
print("5. kesit Qx5 (m3/s): ", Qx5)
Qx6=6*Q52/8/1000
print("6. kesit Qx6 (m3/s): ", Qx6)
Qx7=7*Q52/8/1000
print("7. kesit Qx7 (m3/s): ", Qx7)
Qx8=8*Q52/8/1000
print("8. kesit Qx8 (m3/s): ", Qx8)

rx2 = symbols('rx2')
equation = Eq(A * math.pi * rx2**2 - Qx2 * rx2 - Qx2 * (R3 + t), 0)
rx2_solutions = solve(equation, rx2)
# Pozitif cozumu bulun
positive_rx2 = [sol.evalf() for sol in rx2_solutions if sol > 0]
rx2=round(positive_rx2[0], 4)
print("rx2 (m): ", rx2)
d2=rx2*2*1000
print("d2 (mm): ", d2)

rx3 = symbols('rx3')
equation = Eq(A * math.pi * rx3**2 - Qx3 * rx3 - Qx3 * (R3 + t), 0)
rx3_solutions = solve(equation, rx3)
# Pozitif cozumu bulun
positive_rx3 = [sol.evalf() for sol in rx3_solutions if sol > 0]
rx3=round(positive_rx3[0], 4)
print("rx3 (m): ", rx3)
d3=rx3*2*1000
print("d3 (mm): ", d3)

rx4 = symbols('rx4')
equation = Eq(A * math.pi * rx4**2 - Qx4 * rx4 - Qx4 * (R3 + t), 0)
rx4_solutions = solve(equation, rx4)
# Pozitif cozumu bulun
positive_rx4 = [sol.evalf() for sol in rx4_solutions if sol > 0]
rx4=round(positive_rx4[0], 4)
print("rx4 (m): ", rx4)
d4=rx4*2*1000
print("d4 (mm): ", d4)

rx5 = symbols('rx5')
equation = Eq(A * math.pi * rx5**2 - Qx5 * rx5 - Qx5 * (R3 + t), 0)
rx5_solutions = solve(equation, rx5)
# Pozitif cozumu bulun
positive_rx5 = [sol.evalf() for sol in rx5_solutions if sol > 0]
rx5=round(positive_rx5[0], 4)
print("rx5 (m): ", rx5)
d5=rx5*2*1000
print("d5 (mm): ", d5)

rx6 = symbols('rx6')
equation = Eq(A * math.pi * rx6**2 - Qx6 * rx6 - Qx6 * (R3 + t), 0)
rx6_solutions = solve(equation, rx6)
# Pozitif cozumu bulun
positive_rx6 = [sol.evalf() for sol in rx6_solutions if sol > 0]
rx6=round(positive_rx6[0], 4)
print("rx6 (m): ", rx6)
d6=rx6*2*1000
print("d6 (mm): ", d6)

rx7 = symbols('rx7')
equation = Eq(A * math.pi * rx7**2 - Qx7 * rx7 - Qx7 * (R3 + t), 0)
rx7_solutions = solve(equation, rx7)
# Pozitif cozumu bulun
positive_rx7 = [sol.evalf() for sol in rx7_solutions if sol > 0]
rx7=round(positive_rx7[0], 4)
print("rx7 (m): ", rx7)
d7=rx7*2*1000
print("d7 (mm): ", d7)

rx8 = symbols('rx8')
equation = Eq(A * math.pi * rx8**2 - Qx8 * rx8 - Qx8 * (R3 + t), 0)
rx8_solutions = solve(equation, rx8)
# Pozitif cozumu bulun
positive_rx8 = [sol.evalf() for sol in rx8_solutions if sol > 0]
rx8=round(positive_rx8[0], 4)
print("rx8 (m): ", rx8)
d8=rx8*2*1000
print("d8 (mm): ", d8)

C= round(Qx8/(math.pi*(rx8**2)), 2)
print("Son kesitte su hizi (m/s), C: ", C)

derece = 5
radyan5 = math.radians(derece)
lb=round((200-104.4)/(2*math.tan(radyan5)), 0)
print("lb (mm): ", lb )

print("salyangoz cikis flans capi 150 mm secilirse;" )
db=0.15
Vb=round((4*QT*0.001)/(math.pi*db**2), 2)
print("hiz bulunur (m/s), Vb: ", Vb)

derece2 = 5
radyan52 = math.radians(derece2)
lbI=round((150-104.4)/(2*math.tan(radyan52)), 0)
print("lb' (mm): ", lbI )


