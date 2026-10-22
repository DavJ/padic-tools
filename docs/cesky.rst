Český průvodce
==============

Knihovna už obsahuje numerickou i symbolickou vrstvu a skutečně implementované
propojení s komplexními čísly, kvaterniony, bikvaterniony a geometrickou
algebrou. Je to první použitelná výzkumná verze 0.1; rozsah a další plán
najdeš v :doc:`přehledu vývoje <roadmap>`.

Numerické počítání
------------------

P-adické číslo reprezentujeme jako množinu možných hodnot
``a + O(p^N)``. Známe tedy číslice do určitého řádu, nikoliv nekonečný rozvoj.
Přesnost se přenáší při všech operacích. Dělení může část číslic spotřebovat.
``Qp(5, 20)`` nastavuje výchozí počet významových číslic, nikoliv vždy stejný
absolutní řád chyby.

.. doctest::

   >>> K = Qp(5, 20)
   >>> x = K(Fraction(1, 3))
   >>> (3*x).is_congruent(1)
   True
   >>> sqrt(K(9))**2 == K(9)
   True
   >>> log(exp(K(5))).is_congruent(5)
   True

Při výpočtu ``log`` nebo ``exp`` se kontroluje oblast konvergence. Používá se
Iwasawův logaritmus s výchozí volbou ``log(p)=0``; lze zadat jinou větev.
Pro p=2 platí jiné mezní podmínky než pro lichá prvočísla.

Theta a zeta
------------

``theta3(q, z)`` počítá součet přes všechna celá n ve tvaru
``q^(n*n) * z^n``, tedy multiplikativní theta funkci. Pro ``z=1`` dostaneme
``1 + 2*q + 2*q^4 + ...``. Vyžaduje ``|q|_p < 1``. Zápis
``q=exp(i*pi*tau)`` se mezi komplexním a p-adickým světem automaticky nepřenáší.

``zeta_special(s)`` dává přesnou racionální hodnotu klasické zety pro celé
``s<=0``. ``zeta_p(s, K, branch=r)`` numericky vyhodnocuje konkrétní sudou
Kubotovu–Leopoldtovu větev, i pro necelé p-adické argumenty v Zp. Větev
``r=0`` u p=5 interpoluje hodnoty ``s=1-n`` pro n dělitelné čtyřmi.
Proto samotné ``s=-1`` bez správné větve nestačí k porovnání.

.. doctest::

   >>> zeta_special(-1)
   Fraction(-1, 12)
   >>> zeta_p(-1, Qp(5), branch=2).is_congruent(Fraction(1, 3))
   True
   >>> zeta_p(-3, Qp(5), branch=0).is_congruent(Fraction(-31, 30))
   True

Vazbu tvoří faktor ``(1-p^(n-1))*zeta(1-n)``. Nuly komplexní zety z tohoto
vztahu nelze automaticky ztotožnit s nulami p-adické větve. To zůstává
samostatný matematický problém, nikoliv hotová vlastnost tohoto API.

Transformace a Taylorův rozvoj
------------------------------

``AtomicMeasure`` je konečná vážená míra na Zp s přesnými racionálními
body a vahami. Obsahuje přesné momenty, p-adickou Mellinovu transformaci
na jednotkách a lokální Laplaceovu transformaci s jádrem ``exp(-t*x)``.
Mellinova transformace má explicitní větev kořene jednotky; u p=2
se zvlášť ošetřuje znaménko modulo 4. Omezení míry na jednotky se provádí
výslovně pomocí ``restrict_units()``.

Amiceova transformace používá jádro ``(1+T)^x``. Její binomické momenty
jsou přesné; numerické vyčíslení přidává mez celého zbytku řady.
Na společném oboru platí ``Laplace(t)=Amice(exp(-t)-1)``.
``RationalAmiceMeasure`` navíc zadává míru úplnou racionální Amiceovou
transformací ``P(T)/Q(T)``. Konstruktor ověřuje podmínku omezenosti
koeficientů. Tato vrstva počítá přesné momenty, integrály polynomů,
Laplaceovu transformaci a konvoluci bez konečného seznamu bodů.
Obsahuje také přesné omezení na jednotky a Mellinovu transformaci pro necelé
p-adické exponenty. Ta používá doložený odhad celého Mahlerova zbytku,
zvlášť pro lichá prvočísla a p=2. Při integraci jiného, uživatelského prefixu
Mahlerova rozvoje musí volající doložit mez vynechaných koeficientů.
Inverzní transformace zatím nejsou implementovány.

.. doctest::

   >>> K = Qp(5, 18)
   >>> mira = AtomicMeasure(K, [(1, 1), (2, 2)])
   >>> mira.mellin(6, branch=2).is_congruent(mira.moment(6))
   True
   >>> mira.amice(exp(K(-5))-1).is_congruent(mira.laplace(K(5)))
   True
   >>> rozvoj = taylor(exp, K(5), order=24)
   >>> rozvoj(K(10)).is_congruent(exp(K(10)), precision=18)
   True
   >>> racionalni = RationalAmiceMeasure(K, [1], [2, 1])
   >>> racionalni.moment(3)
   Fraction(1, 8)
   >>> all(racionalni.moment(k) == (-1)**k*(1-2**(k+1))*zeta_special(-k) for k in range(10))
   True
   >>> jednotky = racionalni.restrict_units()
   >>> jednotky.moment(3)
   Fraction(-31, 2)
   >>> jednotky.mellin(3, branch=3).contains(jednotky.moment(3))
   True
   >>> jednotky.mellin(Fraction(1, 2), branch=3).precision
   18

``taylor`` počítá rozvoj ve zvoleném středu. Pro polynomy nad Q dává
přesné racionální koeficienty. Pro p-adický střed podporuje polynomy,
exp, log, log1p, sin a cos; vyčíslení vrací kouli obsahující chybu
koeficientů, vstupu i celého vynechaného zbytku. Počet Taylorových členů
se liší od počtu známých p-adických číslic. Podrobnosti najdeš v
:doc:`kapitole o transformacích <transforms>` a
:doc:`kapitole o Taylorových rozvojích <taylor>`, další příklady pak v
:doc:`kapitole o omezených mírách <measures>` a
:doc:`kapitole o jejich Mellinově transformaci <mellin>`.

Harmonická analýza na Qp
------------------------

``additive_phase`` vrací přesnou fázi p-adického aditivního charakteru.
Jeho komplexní hodnota ``additive_character`` je uložena jako přesný kořen
jednotky; ``.complex(dps=...)`` jej samostatně numericky vyčíslí.
Například p-adická zlomková část čísla 1/6 je pro p=3 rovna 2/3.

``PadicTestFunction`` zadává celou lokálně konstantní funkci s kompaktním
nosičem pomocí hodnot na konečném počtu tříd. Počítá přesný Haarův integrál,
Fourierovu transformaci, inverzi, skalární součin a konvoluci. Používá reálnou
Haarovu míru s objemem Zp rovným jedné. Tato integrace má jiný obor hodnot
než výše popsané p-adické míry. Nosič a mřížka konstantnosti se při Fourierově
transformaci přesně zamění za duální mřížky; knihovna uchovává i tyto údaje.

.. doctest::

   >>> K = Qp(3, 18)
   >>> additive_phase(Fraction(1, 6), K)
   Fraction(2, 3)
   >>> funkce = PadicTestFunction(K, [1, 2, 4])
   >>> funkce.integral().as_fraction()
   Fraction(7, 3)
   >>> spektrum = funkce.fourier()
   >>> spektrum.fourier(inverse=True) == funkce
   True
   >>> spektrum.norm_squared() == funkce.norm_squared()
   True

Zdrojem hodnot jsou racionální čísla a přesné kombinace p-mocninových kořenů
jednotky. U této třídy funkcí jde o konečný přesný výpočet integrálu.
Omezení počtu tříd a řádu kořene vyvolá výslovnou chybu; detaily, znaménka,
normalizaci i odvození najdeš v :doc:`kapitole o harmonické analýze <harmonic>`.

Vlnky, Vladimirovův operátor a difuze
-------------------------------------

``KozyrevBasis`` poskytuje přesný ortogonální rozklad na konečné p-adické
mřížce a jeho rekonstrukci. Vlnky mají amplitudu jedna; přesné normy jsou
uložené zvlášť, takže normalizace nevyžaduje další odmocniny v oboru hodnot.

``vladimirov(f, order=...)`` počítá derivaci na celém Qp pro kladný celý řád.
U funkce s nenulovým integrálem uchovává i celý radiální ocas mimo původní
nosič. ``BallVladimirov`` naopak integruje skoky uvnitř zvolené koule a jeho
difuze zachovává hmotnost. Znaménka, únikový člen a odvození rozdílu obou
operátorů najdeš v :doc:`kapitole o vlnkách a difuzi <wavelets>`.

.. doctest::

   >>> K = Qp(3, 18)
   >>> baze = KozyrevBasis(K)
   >>> signal = PadicTestFunction(K, [1, 2, 4])
   >>> baze.reconstruct(baze.decompose(signal)) == signal
   True
   >>> vlnka = kozyrev_wavelet(K)
   >>> vladimirov(vlnka).as_test_function() == vlnka*3
   True

Příklad ``python examples/padic_diffusion.py`` porovnává jemnější mřížku
se stejnou funkcí a větší kouli s jinou hranicí. Vypisuje spektrum, chybu
hmotnosti, minima hustoty a pokles L2 normy. Vyčíslení difuze používá
komplexní numeriku s nastavenou pracovní přesností, zatímco báze a matice
operátoru zůstávají přesné.

Nerozvětvená kvadratická rozšíření
--------------------------------------------

``UnramifiedQuadraticField`` ověří monický kvadratický polynom s p-adicky
celými přesnými koeficienty a nerozložitelnou redukcí modulo p. Každý prvek
pak nese přesnost celé ideálové koule; valuace je minimum valuací jeho dvou
souřadnic. Uniformizátorem je p a reziduální pole má p na druhou prvků.
K dispozici jsou inverze, algebraická norma, stopa, konjugace a Frobenius.

.. doctest::

   >>> L = UnramifiedQuadraticField(Qp(3, 12), [1, 0, 1], name="i")
   >>> x = L(1, 2)
   >>> (x*x.inverse()).is_congruent(1)
   True
   >>> x.norm().contains(5), x.trace().contains(2)
   (True, True)
   >>> L.residue_field().order
   9

Polynom x^2+1 je vhodný pro Q3(i), ale konstruktor jej odmítne nad Q5,
kde se rozštěpí, i nad Q2, kde definuje rozvětvené rozšíření. Nerozvětvené
kvadratické rozšíření Q2 lze místo toho zadat polynomem x^2+x+1.
``norm()`` vrací algebraickou normu do Qp; p-adickou absolutní hodnotu
vrací ``absolute_value()``. Definice, odvození přesnosti a další příklady
jsou v :doc:`kapitole o rozšířeních <extensions>`. Rozvětvený případ a nové
adaptéry bridge jsou další samostatné kroky.

Jak funguje bridge
------------------

Společné přesné číslo nad Q, případně Q[i], má samostatný komplexní a
p-adický obraz. Uchováváme původní algebraická data a výslovnou volbu místa.
Nesnažíme se převádět konečný p-adický zbytek na údajnou komplexní aproximaci.

.. doctest::

   >>> g = GaussianBridge(Gaussian(Fraction(1, 3), 2))
   >>> local = g.padic_algebra(Qp(3, 20))
   >>> local.parent.d
   Fraction(-1, 1)
   >>> i5 = GaussianBridge(0, 1).split_padic(Qp(5, 20), root_sign=1)
   >>> (i5*i5).is_congruent(-1)
   True

Pro ``p=1 mod 4`` existují dvě volby obrazu i přímo v Qp. Pro ``p=3 mod 4``
potřebujeme kvadratické rozšíření. Nový ověřený engine je popsaný výše;
stávající metoda ``padic_algebra`` nadále vrací koeficientovou algebru.
U p=2 je rozšíření rozvětvené a současná reprezentace neslibuje plnou
aritmetiku s uniformizátorem. Kvadratická algebra je implementována ve všech
těchto případech; rozlišení pole a rozštěpené algebry je podstatné.

Geometrická algebra
-------------------

Hamiltonovy kvaterniony převádíme izomorfně do sudé části Cl(3,0):
``i=e2*e3``, ``j=e3*e1``, ``k=-e1*e2``. Tato identifikace zachovává součiny.
Stejný kód pracuje s racionálními, p-adickými i gaussovskými koeficienty,
a proto i s bikvaterniony. Obsahuje rovněž vektory, bivektory, pseudoskalár,
projekce podle stupně, geometrický a vnější součin.

.. doctest::

   >>> q = Hamilton(1, 2, 3, 4)
   >>> (q*q).to_clifford() == q.to_clifford()*q.to_clifford()
   True
   >>> QuaternionBridge(q).padic(Qp(5)).to_clifford().to_quaternion().is_congruent(QuaternionBridge(q).padic(Qp(5)))
   True

Geometrická algebra tedy sjednocuje algebraické operace. Rozdíl mezi
komplexní a p-adickou konvergencí tím nezaniká. Podrobné definice, odvození
a odkazy jsou v :doc:`kapitole o propojení <bridges>`,
:doc:`kapitole o funkcích <functions>` a :doc:`seznamu zdrojů <references>`.
