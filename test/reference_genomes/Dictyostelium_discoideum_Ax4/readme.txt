Original file "D. discoideum Chromosomal DNA: 1,2,3,4,5,6,M, and floating contigs" (created 05-13-2009 13:53)
File name: dicty_chromosomal.gz

downloaded from

http://dictybase.org/db/cgi-bin/dictyBase/download/blast_databases.pl

# gff files dowloaded from:
http://dictybase.org/Downloads/



reference publication:
Eichinger L, Pachebat JA, Glockner G, Rajandream MA, Sucgang R, Berriman M, Song J, Olsen R, Szafranski K, Xu Q, Tunggal B, Kummerfeld S, Madera M, Konfortov BA, Rivero F, Bankier AT, Lehmann R, Hamlin N, Davies R, Gaudet P, Fey P, Pilcher K, Chen G, Saunders D, Sodergren E, Davis P, Kerhornou A, Nie X, Hall N, Anjard C, Hemphill L, Bason N, Farbrother P, Desany B, Just E, Morio T, Rost R, Churcher C, Cooper J, Haydock S, van Driessche N, Cronin A, Goodhead I, Muzny D, Mourier T, Pain A, Lu M, Harper D, Lindsay R, Hauser H, James K, Quiles M, Madan Babu M, Saito T, Buchrieser C, Wardroper A, Felder M, Thangavelu M, Johnson D, Knights A, Loulseged H, Mungall K, Oliver K, Price C, Quail MA, Urushihara H, Hernandez J, Rabbinowitsch E, Steffen D, Sanders M, Ma J, Kohara Y, Sharp S, Simmonds M, Spiegler S, Tivey A, Sugano S, White B, Walker D, Woodward J, Winckler T, Tanaka Y, Shaulsky G, Schleicher M, Weinstock G, Rosenthal A, Cox EC, Chisholm RL, Gibbs R, Loomis WF, Platzer M, Kay RR, Williams J, Dear PH, Noegel AA, Barrell B, Kuspa A. (2005) The genome of the social amoeba Dictyostelium discoideum. Nature 435(7038): 43-57.

Then Fred Burdet & co have decided to mask the duplicated region of the chromosome, using

bedtools maskfasta -fi dicty_chromosomal2_20090513.fa -fo dicty_chromosomal2_20090513.masked.fa -bed ~/HostPathX/lists/mask_chr2.bed


the latter file containing

DDB0232429      3016083               3768654

(checked, if opened in notepad the masked block is from line 137503-150047)
