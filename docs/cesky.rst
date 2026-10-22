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
potřebujeme kvadratické rozšíření; API pak pracuje s dvojicí koeficientů.
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
