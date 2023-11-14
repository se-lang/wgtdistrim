version 16.1

clear

set seed 42

set obs 1000

generate double w = runiformint(10,50)

replace w =  1 in 1
replace w = 500 in 1000

wgtdistrim w , generate(double wt) upper(.01)

summarize w
scalar mean = r(mean)
scalar sum  = r(sum)

summarize wt

assert reldif(r(mean),mean) < 1e-15
assert reldif(r(sum),sum) < 1e-15

assert r(min) == +1.039a98704a285X+000
assert r(max) == +1.1e76b031c4c33X+006
assert r(sd)  == +1.7e3b0f05ee856X+003

wgtdistrim w , generate(double wt2) upper(.01) lower(0)

assert "`: type wt2'" == "double"
assert wt == wt2

wgtdistrim w , generate(float wt3) upper(.01)

summarize wt3

assert reldif(r(mean),mean) < 1e-8
assert reldif(r(sum),sum) < 1e-8

assert "`: type wt3'" == "float"
assert wt3 != wt // in all obs.

wgtdistrim w , generate(wt4) upper(.01) lower(.01)

assert r(min) == +1.039a980000000X+000
assert r(max) == +1.1e76b00000000X+006
assert r(sd)  == +1.7e3b0f18f1198X+003


capture noisily assert wgtdistrim w , generate(wt) upper(.01)
assert _rc == 111

capture noisily wgtdistrim w , generate(wt5) upper(-.01)
assert _rc == 125

capture noisily wgtdistrim w , generate(wt5) upper(1)
assert _rc == 125

capture noisily wgtdistrim w , generate(wt5) upper(.01) lower(-.01)
assert _rc == 125

capture noisily wgtdistrim w , generate(wt5) upper(.01) lower(1)
assert _rc == 125

capture noisily wgtdistrim w , generate(wt5) upper(.01) iterate(0)
assert _rc == 125

capture noisily wgtdistrim w , generate(wt5) upper(.01) tolerance(-1)
assert _rc == 125

capture noisily wgtdistrim w in 1 , generate(wt5) upper(.01) 
assert _rc == 2001

capture noisily wgtdistrim w if w<1 , generate(wt5) upper(.01) 
assert _rc == 2000

capture noisily wgtdistrim w in 1 , generate(wt5) upper(.01) 
assert _rc == 2001


replace w = 0 in 2
capture noisily wgtdistrim w , generate(wt5) upper(.01)
assert _rc == 459

replace w = -.01 in 2
capture noisily wgtdistrim w , generate(wt5) upper(.01)
assert _rc == 402

replace w = .42 in 2
wgtdistrim w , generate(wt5) upper(.01)


wgtdistrim w in 3/-2 , generate(double wt6) upper(.01)

assert wt6 == w in 3/-2


wgtdistrim w , generate(double w_norm) upper(.01) normalize

summarize w_norm

assert r(mean) == +1.0000000000000X+000
assert r(sd) == +1.8cc2caa43b241X-002
assert r(min) == +1.c37b0d7e366bbX-007
assert r(max) == +1.290fdc11fbaa3X+001
