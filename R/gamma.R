
# let's look at our current 2-parameter gamma distribution for omega
x <- seq(0, 5, length.out=50)
y <- dgamma(x, shape=1.5, rate=3)

plot(x, y, type='l')

y2 <- dgamma(x, shape=2, rate=3)
lines(x, y2, col='red')

# what's the proportion of sites that would have a substantially high omega?
1-pgamma(q=1.5, shape=1.5, rate=3)  # 0.029

1-pgamma(q=1.5, shape=2, rate=3)  # 0.061

# about 20% of sites have w > 1 for shape = 2, rate=3 
# roughly twice as much as before (shape=1.5)
# I suggest bumping shape up to 2 for the LHS runs


# to get the probability of each bin (p vector in INDELIBLE control file)
# INDELIBLE expects the length of this vector to be 1 less
# than the vector of omegas
bins <- diff(pgamma(x, shape=2, rate=3))

# shift omega values over by 1/2 value to get midpoints (w vector in INDELIBLE)
mids <- x + x[2]/2


# Values used in batch_indelible.py
OMEGAS <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65,
          1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.05, 3.15, 3.25, 3.35,
          3.45, 3.55, 3.65, 3.75, 3.85, 3.95, 4.05, 4.15, 4.25, 4.35, 4.45, 4.55, 4.65, 4.75, 4.85, 4.95, 5.05,
          5.15, 5.25, 5.35, 5.45, 5.55, 5.65, 5.75, 5.85, 5.95, 6.05)

# Probabilities of each discrete gamma distro category in OMEGAS
PROP <- c(0.1038610, 0.1432928, 0.1381003, 0.1212857, 0.1020363, 0.0835798, 0.0673901, 0.0535906,
        0.0422005, 0.0329969, 0.0258732, 0.0200207, 0.0154661, 0.0118681, 0.0090903, 0.0070075,
        0.0053782, 0.0040914, 0.0031212, 0.0023785, 0.0017896, 0.0013684, 0.0010189, 0.0007866,
        0.0005856, 0.0004496, 0.0003366, 0.0002510, 0.0001857, 0.0001455, 0.0001097, 0.0000839,
        0.0000653, 0.0000442, 0.0000391, 0.0000294, 0.0000200, 0.0000160, 0.0000124, 0.0000084,
        0.0000066, 0.0000052, 0.0000031, 0.0000022, 0.0000017, 0.0000016, 0.0000009, 0.0000008,
        0.0000008, 0.0000003, 0.0000002, 0.0000002, 0.0000002, 0.0000001, 0.0000003, 0.0000001,
        0.0000002, 0.0000001, 0.0000001, 0.0000001, 0.0000001)


# Attempt to repro those values
repro_omega <- seq(0.0, 6.1, 0.1)
repro_d <- dgamma(repro_omega, shape=1.5, rate=3)
plot(repro_omega, repro_d, type="l")

repro_p <- pgamma(repro_omega, shape=1.5, rate=3)
plot(repro_omega, repro_p)

bin_p <- diff(repro_p)
print(bin_p)

bin_label <- as.factor(repro_omega[1: length(repro_omega)-1] + 0.05)  # the midpoint is the label

plot(bin_label, bin_p)



# Now use the alternate omega distribution where we have more dn/ds that is positively selected
alt_omega <- seq(0.0, 6.1, 0.1)
alt_d <- dgamma(alt_omega, shape=2, rate=3)
plot(alt_omega, alt_d, type="l")

alt_p <- pgamma(alt_omega, shape=2, rate=3)
plot(alt_omega, alt_p)

alt_bin_p <- diff(alt_p)
print(alt_bin_p)

alt_bin_label <- as.factor(alt_omega[1: length(alt_omega)-1] + 0.05)  # the midpoint is the label

plot(alt_bin_label, alt_bin_p)
lines(bin_label, bin_p, col="red")
