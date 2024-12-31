

# Overview ----------------------------------------------------------------


# Sustainable Tapping Model

# David Moore
# December 2024
# davidblakneymoore@gmail.com

# This model is for determining if tapping is sustainable (if annual tree
# growth rates are sufficient for annual tapping). This model makes several
# assumptions, First,the tree is perfectly circular in cross section (and, by
# extension, it is perfectly circular in cross section at all heights, with
# each cross section being perfectly centered on each other and with no vessel
# meandering). Second, each annual increase in the tree's radius is less than
# or equal to the length of the cylindrical portion of the tap hole (which does
# not include the most proximal part of the tap hole that is tapered and that
# was created by the tip of the drill bit). Third, all tap holes are made at
# the same exact height on the tree. Fourth, the heartwood depth is always
# greater than the tap hole depth where the tap hole is installed. Fifth, tap
# holes are made in a perfectly radial direction and parallel with the ground
# and they are directed perfectly to the center of the tree. Sixth, tap holes
# are made in a purely radial direction (there is no axial component to the tap
# hole direction). Seventh, the drill bit's point is not pointy to the extent
# that the cylindrical portion of the drill bit won't completely enter the
# sapwood. Eighth, the cylindrical portion of the tap hole goes deeper than the
# barrel of the tap (the portion of the tap that goes into the tap hole) does
# into the sapwood.


# Packages ----------------------------------------------------------------


if (!require (PBSmapping)) {
  install.packages("PBSmapping")
}
library (PBSmapping)
if (!require (pracma)) {
  install.packages("pracma")
}
library (pracma)
if (!require (sp)) {
  install.packages("sp")
}
library (sp)


# Constants ---------------------------------------------------------------


Minimum_Tapping_Tree_Diameter <- 25 # cm
# Tapping guidelines may suggest minimum diameters of either 10 in or 12 in for
# tapping.
Maximum_Diameter <- 91 # cm
# The maximum diameter at breast height sugar maples can attain is about 36 in
# (https://www.srs.fs.usda.gov/pubs/misc/ag_654/volume_2/acer/saccharum.htm).
Tap_Hole_Depth <- 5.08 # cm
# I'm assuming a 2-in tap hole depth. The 'Tap_Hole_Depth' argument must be
# less than half of the 'Minimum_Tapping_Tree_Diameter' argument.
Tap_Hole_Width <- 0.79375 # cm
# I'm assuming a 0.3125-in tap hole diameter.
Tap_Hole_Drill_Bit_Point_Angle <- (pi / 2) # rad
# Drill bits designed for tapping may have 90 ˚ tips
# (https://www.amazon.com/Liberty-Supply-Professional-Tapping-16/dp/B00IFQKQJQ),
# but they can vary between 90 ˚ and 160 ˚
# (https://www.drillbitwarehouse.com/exploring-the-different-types-of-twist-drill-point-angles/?srsltid=AfmBOoq_KOR8mjd872fkO-GFMnZVxflHdEhg4lf6jpx7B18Cje9Dvkhm).
# The 'Tap_Hole_Drill_Bit_Point_Angle' argument must be between 0 and pi
# (exclusive).
Nonconductive_Wood_Boundary_Width <- 0.5 # cm
Minimum_Tap_Hole_Spacing <- 2 # cm
# This 'Minimum_Tap_Hole_Spacing' argument determines the minimum distance
# between columns of nonconductive wood. I included this parameter in part
# because the point-in-polygon algorithm I use sometimes returns a false
# positive due to rounding.
Tap_Seating_Depth <- 2.54 # cm
# The 'Tap_Seating_Depth' argument specifies how deeply a tap seats in a tap
# hole. Sap is harvested from tissue that is exposed inside the tap hole; it is
# not harvested from tissue that is in contact with the tap.
Tap_Hole_Color <- "black"
Tap_Hole_Wound_Color <- "brown"
Initial_Angle <- 0 # rad
Can_Nonconductive_Wood_Boundaries_Overlap <- TRUE
# The 'Can_Nonconductive_Wood_Boundaries_Overlap' argument specifies if it is
# okay for regions of compartmentalized wood to overlap. Whereas it is always
# bad if an actual tap hole penetrates a nonconductive wood column (whether it
# be compartmentalized tissue or a previous tap hole), it is not known if it is
# bad for the region of compartmentalized wood around a tap hole to overlap (or
# fuse with) another nonconductive wood column.
Distance_Between_Fitted_Points <- 0.01
# The 'Distance_Between_Fitted_Points' argument must be less that half of the
# 'Tap_Hole_Width' argument in order for the function to work (it is best if
# it's much smaller anyway). This argument also must be a positive number.
# Smaller values of this argument lead to more accurate models but also
# increase computational intensity.
Scale_Bar_Text_Shift_Parameter_1 <- 4
# The 'Scale_Bar_Text_Shift_Parameter_1' argument specifies how much the scale
# bar text should shift away from the scale bar when the figure is being stored
# on the hard drive.
Scale_Bar_Text_Shift_Parameter_2 <- 6
# The 'Scale_Bar_Text_Shift_Parameter_2' argument specifies how much the scale
# bar text should shift away from the scale bar when the figure is being made
# in the plotting window in RStudio.
Variable_to_Use_for_Modeling_Growth <- "Annual Radial Growth"
# The 'Variable_to_Use_for_Modeling_Growth' argument must be either "Annual
# Basal Area Increment" or "Annual Radial Growth".
Where_the_Figure_is_Being_Made <- "RStudio Plotting Window"
# The 'Where_the_Figure_is_Being_Made' argument specifies where the figure is
# being made. The two options for this argument are '"RStudio Plotting Window"'
# and '"Hard Drive"'.
Working_Directory <- "/Users/davidblakneymoore/Downloads"
# The 'Working_Directory' argument specifies where the figure is being stored
# on the hard drive (if it's being stored on the hard drive).
Type_of_Distribution_to_Use_for_Modeling_Growth <- "Lognormal"
# The 'Type_of_Distribution_to_Use_for_Modeling_Grouwth' argument must be
# either '"Lognormal"' or '"Uniform"'.
# Additionally, the cylindrical portion of the tap must extend beyond the
# distal portion of the tap barrel into the sapwood.
if ((Tap_Hole_Depth - Tap_Seating_Depth) < (0.5 * Tap_Hole_Width / tan(Tap_Hole_Drill_Bit_Point_Angle / 2))) {
  stop("The cylindrical portion of the tap hole does not extend deeper into the sapwood than the tap barrel.")
}


# Modeling Growth ---------------------------------------------------------


# This section can be refined. It's for determining how to model growth. It's
# based on empirical data about annual radial growth rates and annual basal
# area increments.
if (Variable_to_Use_for_Modeling_Growth == "Annual Radial Growth") {
  if (Type_of_Distribution_to_Use_for_Modeling_Growth == "Lognormal") {
    Lognormal_Distribution_Radial_Increment_Standard_Deviation <- (log(0.3175) - log(0.15875)) / (qnorm(0.75) - qnorm(0.17))
    Lognormal_Distribution_Mean_Radial_Increment = unique(log(0.3175) - (Lognormal_Distribution_Radial_Increment_Standard_Deviation * qnorm(0.75)), log(0.15875) - (Lognormal_Distribution_Radial_Increment_Standard_Deviation * qnorm(0.17)))
  }
} else if (Variable_to_Use_for_Modeling_Growth == "Annual Basal Area Increment") {
  if (Type_of_Distribution_to_Use_for_Modeling_Growth == "Lognormal") {
    Lognormal_Distribution_Basal_Area_Increment_Standard_Deviation <- (log(22.5806) - log(11.61288)) / (qnorm(0.975) - qnorm(0.025))
    Lognormal_Distribution_Mean_Basal_Area_Increment = unique(log(22.5806) - (Intermediary_Lognormal_Distribution_Basal_Area_Increment_Standard_Deviation * qnorm(0.975)), log(11.61288) - (Intermediary_Lognormal_Distribution_Basal_Area_Increment_Standard_Deviation * qnorm(0.025)))
  }
}
# Here, I assumed that 17 % of sugar maple tree ring widths are less than
# 0.15875 cm (0.0625 in), and that 25 % of sugar maple tree ring widths are
# greater than 0.3175 cm (0.125 in)
# (https://mapleresearch.org/wp-content/uploads/wilmot_diameter.pdf). This
# assumption isn't quite accurate because the cited values are for stands and
# not for individual trees. I also assumed that sugar maple tree ring width
# follows a lognormal distribution. My calculations are below.
# Growth could also be modeled with basal area increment following the roughly
# the same protocol. Sugar maple basal area increments are known to vary from
# about 11.61288 cm ^ 2 (1.8 in ^ 2) to about 22.5806 cm ^ 2 (3.5 in ^ 2) for
# healthy dominant or codominant trees
# (https://mapleresearch.org/wp-content/uploads/highyieldgrowth.pdf). I can use
# a uniform distribution to model basal area increments using 11.61288 cm ^ 2
# and 22.5806 cm ^ 2 as the bounds, or I could use a lognormal distribution to
# model the basal area increments if I assume that 95 % of the basal area
# increments are between 11.61288 cm ^ 2 and 22.5806 cm ^ 2.


# The Algorithm -----------------------------------------------------------


Calculating_Areas_of_Overlap_Function <- function (Polygon_of_Interest, Other_Polygons) {
  names(Other_Polygons) <- paste0("Polygon", seq_len(length(Other_Polygons)))
  Modified_Polygon_of_Interest <- data.frame(PID = rep(1, nrow(Polygon_of_Interest)), POS = seq_len(nrow(Polygon_of_Interest)), X = Polygon_of_Interest$Horizontal_Coordinates, Y = Polygon_of_Interest$Vertical_Coordinates)
  Modified_Other_Polygons <- setNames(lapply(seq_len(length(Other_Polygons)), function (x) {
    data.frame(PID = rep((x + 1), nrow(Other_Polygons[[x]])), POS = seq_len(nrow(Other_Polygons[[x]])), X = Other_Polygons[[x]]$Horizontal_Coordinates, Y = Other_Polygons[[x]]$Vertical_Coordinates)
  }), names(Other_Polygons))
  Multipliers <- (-1) ^ ((seq_len(length(Modified_Other_Polygons))) - 1)
  Overlapping_Polygons <- list(NULL)
  Overlapping_Areas <- NULL
  i <- 1
  Flag <- FALSE
  while (Flag == FALSE) {
    if (i == 1) {
      Overlapping_Polygons[[i]] <- lapply(Modified_Other_Polygons, function (x) {
        PBSmapping::joinPolys(Modified_Polygon_of_Interest, x)
      })
      Overlapping_Polygons[[i]] <- Overlapping_Polygons[[i]][which(!sapply(Overlapping_Polygons[[i]], function (x) {
        is.null(x)
      }))]
      if (length(Overlapping_Polygons[[i]]) == 0) {
        Flag <- TRUE
        Overlapping_Area <- 0
      } else if (length(Overlapping_Polygons[[i]]) > 0) {
        Overlapping_Areas[i] <- sum(sapply(Overlapping_Polygons[[i]], function (x) {
          pracma::polyarea(x$X, x$Y)
        }))
      }
    } else if (i == 2) {
      New_Intermediary_Step <- setNames(unlist(lapply(seq_len((length(Modified_Other_Polygons) - 1)), function (k) {
        lapply((k + 1):length(Modified_Other_Polygons), function (j) {
          PBSmapping::joinPolys(Modified_Other_Polygons[[k]], Modified_Other_Polygons[[j]])
        })
      }), recursive = F), unlist(lapply(seq_len((length(Modified_Other_Polygons) - 1)), function (k) {
        lapply((k + 1):length(Modified_Other_Polygons), function (j) {
          paste(names(Modified_Other_Polygons[k]), names(Modified_Other_Polygons[j]), sep = "_")
        })
      }), recursive = F))
      New_Intermediary_Step <- New_Intermediary_Step[which(!(duplicated(lapply(strsplit(names(New_Intermediary_Step), "_"), function (x) {
        y <- as.numeric(gsub("Polygon", "", x))
        y <- y[order(y)]
      }))))]
      New_Intermediary_Step <- New_Intermediary_Step[which(!sapply(New_Intermediary_Step, is.null))]
      New_Intermediary_Step <- lapply(New_Intermediary_Step, function (x) {
        x[((nrow(x)) + 1), ] <- x[1, ]
        x$POS[nrow(x)] <- nrow(x)
        x
      })
      Overlapping_Polygons[[i]] <- lapply(New_Intermediary_Step, function (x) {
        PBSmapping::joinPolys(Modified_Polygon_of_Interest, x)
      })
      Overlapping_Polygons[[i]] <- Overlapping_Polygons[[i]][which(!sapply(Overlapping_Polygons[[i]], is.null))]
      if (length(Overlapping_Polygons[[i]]) == 0) {
        Flag <- TRUE
        Overlapping_Area <- sum(Overlapping_Areas * Multipliers[seq_len(i - 1)])
      } else if (length(Overlapping_Polygons[[i]]) > 0) {
        Overlapping_Areas[i] <- sum(sapply(Overlapping_Polygons[[i]], function (x) {
          pracma::polyarea(x$X, x$Y)
        }))
        Intermediary_Step <- New_Intermediary_Step
      }
    } else if (i > 2) {
      Names_of_Other_Polygons_to_Join <- lapply(lapply(strsplit(names(Intermediary_Step), "_"), function (x) {
        which(!(names(Modified_Other_Polygons) %in% x))
      }), function (y) {
        names(Modified_Other_Polygons)[y]
      })
      New_Intermediary_Step <- unlist(lapply(seq_len(length(Intermediary_Step)), function (k) {
        setNames(lapply(lapply(Names_of_Other_Polygons_to_Join[[k]], function (z) {
          which(names(Modified_Other_Polygons) == z)
        }), function (z) {
          PBSmapping::joinPolys(Intermediary_Step[[k]], Modified_Other_Polygons[[z]])
        }), paste(names(Intermediary_Step[k]), unlist(Names_of_Other_Polygons_to_Join[k]), sep = "_"))
      }), recursive = F)
      New_Intermediary_Step <- New_Intermediary_Step[which(!(duplicated(lapply(strsplit(names(New_Intermediary_Step), "_"), function (x) {
        y <- as.numeric(gsub("Polygon", "", x))
        y <- y[order(y)]
      }))))]
      New_Intermediary_Step <- New_Intermediary_Step[which(!sapply(New_Intermediary_Step, function (x) {
        is.null(x)
      }))]
      New_Intermediary_Step <- lapply(New_Intermediary_Step, function (x) {
        x[((nrow(x)) + 1), ] <- x[1, ]
        x$POS[nrow(x)] <- nrow(x)
        x
      })
      Overlapping_Polygons[[i]] <- lapply(New_Intermediary_Step, function (x) {
        PBSmapping::joinPolys(Modified_Polygon_of_Interest, x)
      })
      Overlapping_Polygons[[i]] <- Overlapping_Polygons[[i]][which(!sapply(Overlapping_Polygons[[i]], is.null))]
      if (length(Overlapping_Polygons[[i]]) == 0) {
        Flag <- TRUE
        Overlapping_Area <- sum(Overlapping_Areas * Multipliers[seq_len(i - 1)])
      } else if (length(Overlapping_Polygons[[i]]) > 0) {
        Overlapping_Areas[i] <- sum(sapply(Overlapping_Polygons[[i]], function (x) {
          pracma::polyarea(x$X, x$Y)
        }))
        Intermediary_Step <- New_Intermediary_Step
      }
    }
    i <- i + 1
    if (i > length(Modified_Other_Polygons)) {
      Overlapping_Area <- sum(Overlapping_Areas * Multipliers[seq_len(i - 1)])
      Flag <- TRUE
    }
  }
  Overlapping_Area
}
Total_Wound_Width <- Tap_Hole_Width + (2 * Nonconductive_Wood_Boundary_Width)
Drill_Bit_Point_Length <- 0.5 * Tap_Hole_Width / tan(Tap_Hole_Drill_Bit_Point_Angle / 2)
Tap_Hole_Wound_Width <- Tap_Hole_Width + (2 * Nonconductive_Wood_Boundary_Width)
Cross_Sectional_Area_of_the_Exposed_Tissue_Inside_Tap_Holes <- (0.5 * Tap_Hole_Width + Drill_Bit_Point_Length) + ((Tap_Hole_Depth - Drill_Bit_Point_Length - Tap_Seating_Depth) * Tap_Hole_Width)
if (Where_the_Figure_is_Being_Made == "Hard Drive") {
  jpeg(paste0(Working_Directory, "/Sustainable Tapping Model Figure.jpeg"), height = 2250, width = 2000)
  par(mar = c(20, 7.5, 5, 1))
  plot(0, axes = F, xlab = "", ylab = "", xlim = c(-((1.05 * (Maximum_Diameter / 2))), ((1.05 * (Maximum_Diameter / 2)))), ylim = c(-((1.05 * (Maximum_Diameter / 2))), ((1.05 * (Maximum_Diameter / 2)))), type = "n", main = "\nTapping Zone Transverse Section\n(as Viewed From Above)", xpd = T, cex.main = 4)
  segments((-((1.1 * (Maximum_Diameter / 2)))), -10, (-(1.1 * (Maximum_Diameter / 2))), 10, xpd = T)
  segments((-((1.1 * (Maximum_Diameter / 2)))), -10, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), -10, xpd = T)
  segments((-((1.1 * (Maximum_Diameter / 2)))), -5, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), -5, xpd = T)
  segments((-((1.1 * (Maximum_Diameter / 2)))), 0, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), 0, xpd = T)
  segments((-((1.1 * (Maximum_Diameter / 2)))), 10, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), 10, xpd = T)
  text(((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_1), -10, "0 cm", xpd = T, cex = 2.5)
  text(((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_1), -5, "5 cm", xpd = T, cex = 2.5)
  text(((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_1), 0, "10 cm", xpd = T, cex = 2.5)
  text(((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_1), 10, "20 cm", xpd = T, cex = 2.5)
  segments(-10, (-((1.1 * (Maximum_Diameter / 2)))), 10, (-(1.1 * (Maximum_Diameter / 2))), xpd = T)
  segments(-10, (-((1.1 * (Maximum_Diameter / 2)))), -10, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), xpd = T)
  segments(-5, (-((1.1 * (Maximum_Diameter / 2)))), -5, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), xpd = T)
  segments(0, (-((1.1 * (Maximum_Diameter / 2)))), 0, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), xpd = T)
  segments(10, (-((1.1 * (Maximum_Diameter / 2)))), 10, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), xpd = T)
  text(-10, ((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_1), "0 cm", xpd = T, srt = 270, cex = 2.5)
  text(-5, ((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_1), "5 cm", xpd = T, srt = 270, cex = 2.5)
  text(0, ((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_1), "10 cm", xpd = T, srt = 270, cex = 2.5)
  text(10, ((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_1), "20 cm", xpd = T, srt = 270, cex = 2.5)
  legend("bottom", xpd = T, inset = c(0, -0.125), legend = c("The Tap Hole", "Compartmentalized Wood Around the Tap Hole"), col = c(Tap_Hole_Color, Tap_Hole_Wound_Color), pch = 15, cex = 2.5)
} else if (Where_the_Figure_is_Being_Made == "RStudio Plotting Window") {
  par(mar = c(7, 4, 4, 1))
  plot(0, axes = F, xlab = "", ylab = "", xlim = c(-((1.05 * (Maximum_Diameter / 2))), ((1.05 * (Maximum_Diameter / 2)))), ylim = c(-((1.05 * (Maximum_Diameter / 2))), ((1.05 * (Maximum_Diameter / 2)))), type = "n", main = "\nTapping Zone Transverse Section\n(as Viewed From Above)", xpd = T, cex.main = 2)
  segments((-((1.1 * (Maximum_Diameter / 2)))), -10, (-(1.1 * (Maximum_Diameter / 2))), 10, xpd = T)
  segments((-((1.1 * (Maximum_Diameter / 2)))), -10, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), -10, xpd = T)
  segments((-((1.1 * (Maximum_Diameter / 2)))), -5, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), -5, xpd = T)
  segments((-((1.1 * (Maximum_Diameter / 2)))), 0, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), 0, xpd = T)
  segments((-((1.1 * (Maximum_Diameter / 2)))), 10, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), 10, xpd = T)
  text(((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_2), -10, "0 cm", xpd = T)
  text(((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_2), -5, "5 cm", xpd = T)
  text(((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_2), 0, "10 cm", xpd = T)
  text(((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_2), 10, "20 cm", xpd = T)
  segments(-10, (-((1.1 * (Maximum_Diameter / 2)))), 10, (-(1.1 * (Maximum_Diameter / 2))), xpd = T)
  segments(-10, (-((1.1 * (Maximum_Diameter / 2)))), -10, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), xpd = T)
  segments(-5, (-((1.1 * (Maximum_Diameter / 2)))), -5, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), xpd = T)
  segments(0, (-((1.1 * (Maximum_Diameter / 2)))), 0, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), xpd = T)
  segments(10, (-((1.1 * (Maximum_Diameter / 2)))), 10, ((-(1.1 * (Maximum_Diameter / 2))) - 1.25), xpd = T)
  text(-10, ((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_2), "0 cm", xpd = T, srt = 270)
  text(-5, ((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_2), "5 cm", xpd = T, srt = 270)
  text(0, ((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_2), "10 cm", xpd = T, srt = 270)
  text(10, ((-(1.1 * (Maximum_Diameter / 2))) - Scale_Bar_Text_Shift_Parameter_2), "20 cm", xpd = T, srt = 270)
  legend("bottom", xpd = T, inset = c(0, -0.2), legend = c("The Tap Hole", "Compartmentalized Wood Around the Tap Hole"), col = c(Tap_Hole_Color, Tap_Hole_Wound_Color), pch = 15)
}
Radii <- NULL
Angles <- NULL
Has_the_First_Tap_Been_Installed_Yet <- FALSE
Tap_Hole_Polygon_Coordinates <- list(NULL)
Tap_Hole_Wound_Polygon_Coordinates <- list(NULL)
Compartmentalized_Wood_Only <- list(NULL)
Exposed_Tissue_Inside_the_Tap_Hole <- list(NULL)
i <- 1
repeat {
  if (i == 1) {
    if (Variable_to_Use_for_Modeling_Growth == "Annual Radial Growth") {
      if (Type_of_Distribution_to_Use_for_Modeling_Growth == "Lognormal") {
        Radii[i] <- rlnorm(1, Lognormal_Distribution_Mean_Radial_Increment, Lognormal_Distribution_Radial_Increment_Standard_Deviation)
      } else if (Type_of_Distribution_to_Use_for_Modeling_Growth == "Uniform") {
        Radii[i] <- runif(1, 0.15875, 0.3175)
      }
    } else if (Variable_to_Use_for_Modeling_Growth == "Annual Basal Area Increment") {
      if (Type_of_Distribution_to_Use_for_Modeling_Growth == "Lognormal") {
        Radii[i] <- sqrt(rlnorm(1, Lognormal_Distribution_Mean_Basal_Area_Increment, Lognormal_Distribution_Basal_Area_Increment_Standard_Deviation) / pi)
      } else if (Type_of_Distribution_to_Use_for_Modeling_Growth == "Uniform") {
        Radii[i] <- sqrt(runif(1, 11.61288, 22.5806) / pi)
      }
    }
  } else if (i > 1) {
    if (Variable_to_Use_for_Modeling_Growth == "Annual Radial Growth") {
      if (Type_of_Distribution_to_Use_for_Modeling_Growth == "Lognormal") {
        Radii[i] <- Radii[(i - 1)] + rlnorm(1, Lognormal_Distribution_Mean_Radial_Increment, Lognormal_Distribution_Radial_Increment_Standard_Deviation)
      } else if (Type_of_Distribution_to_Use_for_Modeling_Growth == "Uniform") {
        Radii[i] <- Radii[(i - 1)] + runif(1, 0.15875, 0.3175)
      }
    } else if (Variable_to_Use_for_Modeling_Growth == "Annual Basal Area Increment") {
      if (Type_of_Distribution_to_Use_for_Modeling_Growth == "Lognormal") {
        Radii[i] <- sqrt(((pi * (Radii[(i - 1)] ^ 2)) + rlnorm(1, Lognormal_Distribution_Mean_Basal_Area_Increment, Lognormal_Distribution_Basal_Area_Increment_Standard_Deviation)) / pi)
      } else if (Type_of_Distribution_to_Use_for_Modeling_Growth == "Uniform") {
        Radii[i] <- sqrt(((pi * (Radii[(i - 1)] ^ 2)) + runif(1, 11.61288, 22.5806)) / pi)
      }
    }
  }
  if (Radii[i] > (Maximum_Diameter / 2)) {
    # warning ("The maximum sugar maple diameter at breast height has been reached.")
    break
  }
  if (i == 1) {
    Data_Frame <- data.frame(Year = i, Diameter = (2 * Radii[i]), Tap_Hole_Made = ifelse((Radii[i] > (Minimum_Tapping_Tree_Diameter / 2)), T, F), Area_of_Overlap_Between_This_Tap_Hole_and_Other_Nonconductive_Wood_Columns = NA, Proportion_of_Exposed_Tissue_in_the_Tap_Hole_That_Overlaps_With_Other_Nonconductive_Wood_Columns = NA)
  } else if (i > 1) {
    Data_Frame[i, ] <- data.frame(Year = i, Diameter = (2 * Radii[i]), Tap_Hole_Made = ifelse((Radii[i] > (Minimum_Tapping_Tree_Diameter / 2)), T, F), Area_of_Overlap_Between_This_Tap_Hole_and_Other_Nonconductive_Wood_Columns = NA, Proportion_of_Exposed_Tissue_in_the_Tap_Hole_That_Overlaps_With_Other_Nonconductive_Wood_Columns = NA)
  }
  Angle_Increment_1 <- Distance_Between_Fitted_Points / Radii[i]
  Tree_Ring_Coordinates <- as.data.frame(t(as.matrix(as.data.frame(lapply(seq_len(floor((2 * pi * Radii[i]) / Distance_Between_Fitted_Points)), function (l) {
    c(Radii[i] * c(cos((l - 1) * Angle_Increment_1), sin((l - 1) * Angle_Increment_1)))
  })))))
  rownames(Tree_Ring_Coordinates) <- NULL
  colnames(Tree_Ring_Coordinates) <- c("Horizontal_Coordinates", "Vertical_Coordinates")
  Tree_Ring_Coordinates[((nrow(Tree_Ring_Coordinates)) + 1), ] <- Tree_Ring_Coordinates[1, ]
  polygon(Tree_Ring_Coordinates$Horizontal_Coordinates, Tree_Ring_Coordinates$Vertical_Coordinates)
  if (Radii[i] > (Minimum_Tapping_Tree_Diameter / 2)) {
    Sagitta_Length <- Radii[i] - sqrt((Radii[i] ^ 2) - ((Tap_Hole_Width ^ 2) / 4))
    if (Has_the_First_Tap_Been_Installed_Yet == FALSE) {
      k <- 1
      Angles[k] <- Initial_Angle
    } else if (Has_the_First_Tap_Been_Installed_Yet == TRUE) {
      k <- k + 1
      Reference_Radius <- Radii[i] - (Tap_Hole_Depth - Drill_Bit_Point_Length)
      Reference_Hypotenuse <- sqrt((Reference_Radius ^ 2) + ((Tap_Hole_Width / 2) ^ 2))
      if (Can_Nonconductive_Wood_Boundaries_Overlap == T) {
        Spacing <- (Tap_Hole_Width / 2) + Nonconductive_Wood_Boundary_Width + Minimum_Tap_Hole_Spacing
      } else if (Can_Nonconductive_Wood_Boundaries_Overlap == F) {
        Spacing <- (Tap_Hole_Width / 2) + (2 * Nonconductive_Wood_Boundary_Width) + Minimum_Tap_Hole_Spacing
      }
      New_Distance <- sqrt((Reference_Hypotenuse ^ 2) - (Spacing ^ 2))
      Angles[k] <- Angles[(k - 1)] + atan(Spacing / New_Distance) + atan((Tap_Hole_Width / 2) / Reference_Radius)
    }
    Point_1 <- (Radii[i] - Sagitta_Length) * c(cos(Angles[k]), sin(Angles[k]))
    Point_2 <- (Radii[i] - Tap_Hole_Depth) * c(cos(Angles[k]), sin(Angles[k]))
    Counterclockwise_Shift <- (Tap_Hole_Width / 2) * c(-sin(Angles[k]), cos(Angles[k]))
    Clockwise_Shift <- (Tap_Hole_Width / 2) * c(sin(Angles[k]), -cos(Angles[k]))
    Distal_Shift <- Drill_Bit_Point_Length * c(cos(Angles[k]), sin(Angles[k]))
    Point_3 <- Point_1 + Counterclockwise_Shift
    Point_4 <- Point_1 + Clockwise_Shift
    Point_5 <- Point_2 + Counterclockwise_Shift + Distal_Shift
    Point_6 <- Point_2 + Clockwise_Shift + Distal_Shift
    Point_7 <- c(Radii[i] * cos(Angles[k]), Radii[i] * sin(Angles[k]))
    if (Has_the_First_Tap_Been_Installed_Yet == TRUE) {
      Tap_Seating_Depth_Radius <- Radii[i] - Tap_Seating_Depth
      Point_8 <- Tap_Seating_Depth_Radius * c(cos(Angles[k]), sin(Angles[k]))
      Point_9 <- Point_8 + Counterclockwise_Shift
      Point_10 <- Point_8 + Clockwise_Shift
      Exposed_Tissue_Inside_the_Tap_Hole[[k]] <- data.frame(Horizontal_Coordinates = c(Point_10[1], Point_6[1], Point_2[1], Point_5[1], Point_9[1], Point_10[1]), Vertical_Coordinates = c(Point_10[2], Point_6[2], Point_2[2], Point_5[2], Point_9[2], Point_10[2]))
    }
    Distance_Between_Points_3_and_4_and_Point_7 <- c(unique(abs(Point_7[1] - Point_3[1]), abs(Point_7[1] - Point_4[1])), unique(abs(Point_7[2] - Point_3[2]), abs(Point_7[2] - Point_4[2])))
    Points_on_the_Perimeter_That_Are_Part_of_the_Tap_Hole <- Tree_Ring_Coordinates[which((abs(Tree_Ring_Coordinates$Horizontal_Coordinates - Point_7[1]) < Distance_Between_Points_3_and_4_and_Point_7[1]) & (abs(Tree_Ring_Coordinates$Vertical_Coordinates - Point_7[2]) < Distance_Between_Points_3_and_4_and_Point_7[2])), ]
    Tap_Hole_Polygon_Coordinates[[k]] <- data.frame(Horizontal_Coordinates = c(Point_4[1], Point_6[1], Point_2[1], Point_5[1], Point_3[1], Points_on_the_Perimeter_That_Are_Part_of_the_Tap_Hole$Horizontal_Coordinates, Point_4[1]), Vertical_Coordinates = c(Point_4[2], Point_6[2], Point_2[2], Point_5[2], Point_3[2], Points_on_the_Perimeter_That_Are_Part_of_the_Tap_Hole$Vertical_Coordinates, Point_4[2]))
    polygon(Tap_Hole_Polygon_Coordinates[[k]]$Horizontal_Coordinates, Tap_Hole_Polygon_Coordinates[[k]]$Vertical_Coordinates, col = Tap_Hole_Color)
    Angle_Increment_2 <- Distance_Between_Fitted_Points / Nonconductive_Wood_Boundary_Width
    Point_2_Circle_Coordinates <- as.data.frame(t(as.matrix(as.data.frame(lapply(seq_len(floor((2 * pi * Nonconductive_Wood_Boundary_Width) / Distance_Between_Fitted_Points)), function (l) {
      c(Nonconductive_Wood_Boundary_Width * c(cos((l - 1) * Angle_Increment_2), sin((l - 1) * Angle_Increment_2))) + Point_2
    })))))
    rownames(Point_2_Circle_Coordinates) <- NULL
    colnames(Point_2_Circle_Coordinates) <- c("Horizontal_Coordinates", "Vertical_Coordinates")
    Point_2_Circle_Coordinates[((nrow(Point_2_Circle_Coordinates)) + 1), ] <- Point_2_Circle_Coordinates[1, ]
    Point_5_Circle_Coordinates <- as.data.frame(t(as.matrix(as.data.frame(lapply(seq_len(floor((2 * pi * Nonconductive_Wood_Boundary_Width) / Distance_Between_Fitted_Points)), function (l) {
      c(Nonconductive_Wood_Boundary_Width * c(cos((l - 1) * Angle_Increment_2), sin((l - 1) * Angle_Increment_2))) + Point_5
    })))))
    rownames(Point_5_Circle_Coordinates) <- NULL
    colnames(Point_5_Circle_Coordinates) <- c("Horizontal_Coordinates", "Vertical_Coordinates")
    Point_5_Circle_Coordinates[((nrow(Point_5_Circle_Coordinates)) + 1), ] <- Point_5_Circle_Coordinates[1, ]
    Point_6_Circle_Coordinates <- as.data.frame(t(as.matrix(as.data.frame(lapply(seq_len(floor((2 * pi * Nonconductive_Wood_Boundary_Width) / Distance_Between_Fitted_Points)), function (l) {
      c(Nonconductive_Wood_Boundary_Width * c(cos((l - 1) * Angle_Increment_2), sin((l - 1) * Angle_Increment_2))) + Point_6
    })))))
    rownames(Point_6_Circle_Coordinates) <- NULL
    colnames(Point_6_Circle_Coordinates) <- c("Horizontal_Coordinates", "Vertical_Coordinates")
    Point_6_Circle_Coordinates[((nrow(Point_6_Circle_Coordinates)) + 1), ] <- Point_6_Circle_Coordinates[1, ]
    Point_5a_Angle <- Angles[k] + (pi / 2)
    Point_5a_Shift <- Nonconductive_Wood_Boundary_Width * c(cos(Point_5a_Angle), sin(Point_5a_Angle))
    Point_5a <- (Point_5 + Point_5a_Shift)
    Point_5b_Angle <- Angles[k] + (pi / 2) + (Tap_Hole_Drill_Bit_Point_Angle / 2)
    Point_5b_Shift <- Nonconductive_Wood_Boundary_Width * c(cos(Point_5b_Angle), sin(Point_5b_Angle))
    Point_5b <- (Point_5 + Point_5b_Shift)
    Point_6a_Angle <- Angles[k] - (pi / 2)
    Point_6a_Shift <- Nonconductive_Wood_Boundary_Width * c(cos(Point_6a_Angle), sin(Point_6a_Angle))
    Point_6a <- (Point_6 + Point_6a_Shift)
    Point_6b_Angle <- Angles[k] - (pi / 2) - (Tap_Hole_Drill_Bit_Point_Angle / 2)
    Point_6b_Shift <- Nonconductive_Wood_Boundary_Width * c(cos(Point_6b_Angle), sin(Point_6b_Angle))
    Point_6b <- (Point_6 + Point_6b_Shift)
    Point_2a_Angle <- Angles[k] + (pi / 2) + (Tap_Hole_Drill_Bit_Point_Angle / 2)
    Point_2a_Shift <- Nonconductive_Wood_Boundary_Width * c(cos(Point_2a_Angle), sin(Point_2a_Angle))
    Point_2a <- (Point_2 + Point_2a_Shift)
    Point_2b_Angle <- Angles[k] - (pi / 2) - (Tap_Hole_Drill_Bit_Point_Angle / 2)
    Point_2b_Shift <- Nonconductive_Wood_Boundary_Width * c(cos(Point_2b_Angle), sin(Point_2b_Angle))
    Point_2b <- (Point_2 + Point_2b_Shift)
    Circle_Centered_at_the_Origin_Coordinates <- as.data.frame(t(as.matrix(as.data.frame(lapply(seq_len(floor((2 * pi * Nonconductive_Wood_Boundary_Width) / Distance_Between_Fitted_Points)), function (l) {
      c(Nonconductive_Wood_Boundary_Width * c(cos((l - 1) * Angle_Increment_2), sin((l - 1) * Angle_Increment_2)))
    })))))
    rownames(Circle_Centered_at_the_Origin_Coordinates) <- NULL
    colnames(Circle_Centered_at_the_Origin_Coordinates) <- c("Horizontal_Coordinates", "Vertical_Coordinates")
    Circle_Centered_at_the_Origin_Coordinates[((nrow(Circle_Centered_at_the_Origin_Coordinates)) + 1), ] <- Circle_Centered_at_the_Origin_Coordinates[1, ]
    if ((cos(Angles[k]) >= 0) & (sin(Angles[k]) >= 0)) {
      Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_2 <- ((Circle_Centered_at_the_Origin_Coordinates$Horizontal_Coordinates < Point_2b_Shift[1]) & (Circle_Centered_at_the_Origin_Coordinates$Vertical_Coordinates < Point_2a_Shift[2]))
    } else if ((cos(Angles[k]) < 0) & (sin(Angles[k]) >= 0)) {
      Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_2 <- ((Circle_Centered_at_the_Origin_Coordinates$Horizontal_Coordinates > Point_2a_Shift[1]) & (Circle_Centered_at_the_Origin_Coordinates$Vertical_Coordinates < Point_2b_Shift[2]))
    } else if ((cos(Angles[k]) < 0) & (sin(Angles[k]) < 0)) {
      Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_2 <- ((Circle_Centered_at_the_Origin_Coordinates$Horizontal_Coordinates > Point_2b_Shift[1]) & (Circle_Centered_at_the_Origin_Coordinates$Vertical_Coordinates > Point_2a_Shift[2]))
    } else if ((cos(Angles[k]) >= 0) & (sin(Angles[k]) < 0)) {
      Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_2 <- ((Circle_Centered_at_the_Origin_Coordinates$Horizontal_Coordinates < Point_2a_Shift[1]) & (Circle_Centered_at_the_Origin_Coordinates$Vertical_Coordinates > Point_2b_Shift[2]))
    }
    Run_Length_Encoding <- rle(Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_2)
    Run_Length_Encoding <- data.frame(Lengths = Run_Length_Encoding$lengths, Values = Run_Length_Encoding$values)
    if (length(which(Run_Length_Encoding$Values == TRUE)) == 1) {
      Point_2_Circle_Coordinates <- Point_2_Circle_Coordinates[which(Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_2 == T), ]
    } else if (length(which(Run_Length_Encoding$Values == TRUE)) == 2) {
      Cumulative_Sums <- cumsum(Run_Length_Encoding$Lengths)
      Run_Length_Encoding$Starting_Values <- c(1, (Cumulative_Sums[seq_len((length(Cumulative_Sums)) - 1)]) + 1)
      Run_Length_Encoding$Ending_Values <- Cumulative_Sums
      Point_2_Circle_Coordinates <- rbind(Point_2_Circle_Coordinates[c(Run_Length_Encoding$Starting_Values[3]:Run_Length_Encoding$Ending_Values[3]), ], Point_2_Circle_Coordinates[c(Run_Length_Encoding$Starting_Values[1]:Run_Length_Encoding$Ending_Values[1]), ])
      rownames(Point_2_Circle_Coordinates) <- NULL
    }
    Point_2_Circle_Coordinates <- Point_2_Circle_Coordinates[which(!duplicated(Point_2_Circle_Coordinates)), ]
    Point_5_Bisecting_Vector <- Point_5a_Shift + Point_5b_Shift
    Point_5_Bisecting_Vector <- Point_5_Bisecting_Vector / sqrt((Point_5_Bisecting_Vector[1] ^ 2) + (Point_5_Bisecting_Vector[2] ^ 2))
    Point_5_Bisecting_Angle <- atan2(Point_5_Bisecting_Vector[2], Point_5_Bisecting_Vector[1])
    if ((cos(Point_5_Bisecting_Angle) >= 0) & (sin(Point_5_Bisecting_Angle) >= 0)) {
      Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_5 <- ((Circle_Centered_at_the_Origin_Coordinates$Horizontal_Coordinates > Point_5b_Shift[1]) & (Circle_Centered_at_the_Origin_Coordinates$Vertical_Coordinates > Point_5a_Shift[2]))
    } else if ((cos(Point_5_Bisecting_Angle) < 0) & (sin(Point_5_Bisecting_Angle) >= 0)) {
      Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_5 <- ((Circle_Centered_at_the_Origin_Coordinates$Horizontal_Coordinates < Point_5a_Shift[1]) & (Circle_Centered_at_the_Origin_Coordinates$Vertical_Coordinates > Point_5b_Shift[2]))
    } else if ((cos(Point_5_Bisecting_Angle) < 0) & (sin(Point_5_Bisecting_Angle) < 0)) {
      Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_5 <- ((Circle_Centered_at_the_Origin_Coordinates$Horizontal_Coordinates < Point_5b_Shift[1]) & (Circle_Centered_at_the_Origin_Coordinates$Vertical_Coordinates < Point_5a_Shift[2]))
    } else if ((cos(Point_5_Bisecting_Angle) >= 0) & (sin(Point_5_Bisecting_Angle) < 0)) {
      Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_5 <- ((Circle_Centered_at_the_Origin_Coordinates$Horizontal_Coordinates > Point_5a_Shift[1]) & (Circle_Centered_at_the_Origin_Coordinates$Vertical_Coordinates < Point_5b_Shift[2]))
    }
    Run_Length_Encoding <- rle(Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_5)
    Run_Length_Encoding <- data.frame(Lengths = Run_Length_Encoding$lengths, Values = Run_Length_Encoding$values)
    if (length(which(Run_Length_Encoding$Values == TRUE)) == 1) {
      Point_5_Circle_Coordinates <- Point_5_Circle_Coordinates[which(Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_5 == T), ]
    } else if (length(which(Run_Length_Encoding$Values == TRUE)) == 2) {
      Cumulative_Sums <- cumsum(Run_Length_Encoding$Lengths)
      Run_Length_Encoding$Starting_Values <- c(1, (Cumulative_Sums[seq_len((length(Cumulative_Sums)) - 1)]) + 1)
      Run_Length_Encoding$Ending_Values <- Cumulative_Sums
      Point_5_Circle_Coordinates <- rbind(Point_5_Circle_Coordinates[c(Run_Length_Encoding$Starting_Values[3]:Run_Length_Encoding$Ending_Values[3]), ], Point_5_Circle_Coordinates[c(Run_Length_Encoding$Starting_Values[1]:Run_Length_Encoding$Ending_Values[1]), ])
      rownames(Point_5_Circle_Coordinates) <- NULL
    }
    Point_5_Circle_Coordinates <- Point_5_Circle_Coordinates[which(!duplicated(Point_5_Circle_Coordinates)), ]
    Point_6_Bisecting_Vector <- Point_6a_Shift + Point_6b_Shift
    Point_6_Bisecting_Vector <- Point_6_Bisecting_Vector / sqrt((Point_6_Bisecting_Vector[1] ^ 2) + (Point_6_Bisecting_Vector[2] ^ 2))
    Point_6_Bisecting_Angle <- atan2(Point_6_Bisecting_Vector[2], Point_6_Bisecting_Vector[1])
    if ((cos(Point_6_Bisecting_Angle) >= 0) & (sin(Point_6_Bisecting_Angle) >= 0)) {
      Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_6 <- ((Circle_Centered_at_the_Origin_Coordinates$Horizontal_Coordinates > Point_6a_Shift[1]) & (Circle_Centered_at_the_Origin_Coordinates$Vertical_Coordinates > Point_6b_Shift[2]))
    } else if ((cos(Point_6_Bisecting_Angle) < 0) & (sin(Point_6_Bisecting_Angle) >= 0)) {
      Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_6 <- ((Circle_Centered_at_the_Origin_Coordinates$Horizontal_Coordinates < Point_6b_Shift[1]) & (Circle_Centered_at_the_Origin_Coordinates$Vertical_Coordinates > Point_6a_Shift[2]))
    } else if ((cos(Point_6_Bisecting_Angle) < 0) & (sin(Point_6_Bisecting_Angle) < 0)) {
      Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_6 <- ((Circle_Centered_at_the_Origin_Coordinates$Horizontal_Coordinates < Point_6a_Shift[1]) & (Circle_Centered_at_the_Origin_Coordinates$Vertical_Coordinates < Point_6b_Shift[2]))
    } else if ((cos(Point_6_Bisecting_Angle) >= 0) & (sin(Point_6_Bisecting_Angle) < 0)) {
      Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_6 <- ((Circle_Centered_at_the_Origin_Coordinates$Horizontal_Coordinates > Point_6b_Shift[1]) & (Circle_Centered_at_the_Origin_Coordinates$Vertical_Coordinates < Point_6a_Shift[2]))
    }
    Run_Length_Encoding <- rle(Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_6)
    Run_Length_Encoding <- data.frame(Lengths = Run_Length_Encoding$lengths, Values = Run_Length_Encoding$values)
    if (length(which(Run_Length_Encoding$Values == TRUE)) == 1) {
      Point_6_Circle_Coordinates <- Point_6_Circle_Coordinates[which(Circle_Centered_at_the_Origin_Coordinates$Are_Points_in_the_Arc_for_Near_Point_6 == T), ]
    } else if (length(which(Run_Length_Encoding$Values == TRUE)) == 2) {
      Cumulative_Sums <- cumsum(Run_Length_Encoding$Lengths)
      Run_Length_Encoding$Starting_Values <- c(1, (Cumulative_Sums[seq_len((length(Cumulative_Sums)) - 1)]) + 1)
      Run_Length_Encoding$Ending_Values <- Cumulative_Sums
      Point_6_Circle_Coordinates <- rbind(Point_6_Circle_Coordinates[c(Run_Length_Encoding$Starting_Values[3]:Run_Length_Encoding$Ending_Values[3]), ], Point_6_Circle_Coordinates[c(Run_Length_Encoding$Starting_Values[1]:Run_Length_Encoding$Ending_Values[1]), ])
      rownames(Point_6_Circle_Coordinates) <- NULL
    }
    Point_6_Circle_Coordinates <- Point_6_Circle_Coordinates[which(!duplicated(Point_6_Circle_Coordinates)), ]
    New_Sagitta_Length <- Radii[i] - sqrt((Radii[i] ^ 2) - ((Tap_Hole_Wound_Width ^ 2) / 4))
    New_Point_1 <- (Radii[i] - New_Sagitta_Length) * c(cos(Angles[k]), sin(Angles[k]))
    Counterclockwise_Shift <- (Tap_Hole_Wound_Width / 2) * c(-sin(Angles[k]), cos(Angles[k]))
    Clockwise_Shift <- (Tap_Hole_Wound_Width / 2) * c(sin(Angles[k]), -cos(Angles[k]))
    New_Point_3 <- New_Point_1 + Counterclockwise_Shift
    New_Point_4 <- New_Point_1 + Clockwise_Shift
    Distance_Between_New_Points_3_and_4_and_Point_7 <- c(unique(abs(Point_7[1] - New_Point_3[1]), abs(Point_7[1] - New_Point_4[1])), unique(abs(Point_7[2] - New_Point_3[2]), abs(Point_7[2] - New_Point_4[2])))
    Points_on_the_Perimeter_That_Are_Part_of_the_Tap_Hole_Wound <- Tree_Ring_Coordinates[which((abs(Tree_Ring_Coordinates$Horizontal_Coordinates - Point_7[1]) < Distance_Between_New_Points_3_and_4_and_Point_7[1]) & (abs(Tree_Ring_Coordinates$Vertical_Coordinates - Point_7[2]) < Distance_Between_New_Points_3_and_4_and_Point_7[2])), ]
    Tap_Hole_Wound_Polygon_Coordinates[[k]] <- data.frame(Horizontal_Coordinates = c(New_Point_4[1], Point_6a[1], rev(Point_6_Circle_Coordinates$Horizontal_Coordinates), Point_6b[1], Point_2b[1], rev(Point_2_Circle_Coordinates$Horizontal_Coordinates), Point_2a[1], Point_5b[1], rev(Point_5_Circle_Coordinates$Horizontal_Coordinates), Point_5a[1], New_Point_3[1], Points_on_the_Perimeter_That_Are_Part_of_the_Tap_Hole_Wound$Horizontal_Coordinates, New_Point_4[1]), Vertical_Coordinates = c(New_Point_4[2], Point_6a[2], rev(Point_6_Circle_Coordinates$Vertical_Coordinates), Point_6b[2], Point_2b[2], rev(Point_2_Circle_Coordinates$Vertical_Coordinates), Point_2a[2], Point_5b[2], rev(Point_5_Circle_Coordinates$Vertical_Coordinates), Point_5a[2], New_Point_3[2], Points_on_the_Perimeter_That_Are_Part_of_the_Tap_Hole_Wound$Vertical_Coordinates, New_Point_4[2]))
    Indices <- which(sapply(seq_len(nrow(Tap_Hole_Wound_Polygon_Coordinates[[k]])), function (x) {
      any(sapply(seq_len(nrow(Tap_Hole_Polygon_Coordinates[[k]])), function (y) {
        all(as.vector(unlist(Tap_Hole_Wound_Polygon_Coordinates[[k]][x, ])) == as.vector(unlist(Tap_Hole_Polygon_Coordinates[[k]][y, ])))
      }))
    }))
    Tap_Hole_Wound_Polygon_Coordinates_First_Subset <- Tap_Hole_Wound_Polygon_Coordinates[[k]][seq_len((min(Indices)) - 1), ]
    Tap_Hole_Wound_Polygon_Coordinates_Last_Subset <- Tap_Hole_Wound_Polygon_Coordinates[[k]][((max(Indices)) + 1):nrow(Tap_Hole_Wound_Polygon_Coordinates[[k]]), ]
    Tap_Hole_Wound_Polygon_Coordinates_Middle_Subset <- rbind(Point_3, Point_5, Point_2, Point_6, Point_4)
    rownames(Tap_Hole_Wound_Polygon_Coordinates_Middle_Subset) <- NULL
    Tap_Hole_Wound_Polygon_Coordinates_Middle_Subset <- as.data.frame(Tap_Hole_Wound_Polygon_Coordinates_Middle_Subset)
    colnames(Tap_Hole_Wound_Polygon_Coordinates_Middle_Subset) <- colnames(Tap_Hole_Wound_Polygon_Coordinates[[k]])
    Compartmentalized_Wood_Only[[k]] <- rbind(Tap_Hole_Wound_Polygon_Coordinates_First_Subset, Tap_Hole_Wound_Polygon_Coordinates_Middle_Subset, Tap_Hole_Wound_Polygon_Coordinates_Last_Subset)
    polygon(Compartmentalized_Wood_Only[[k]]$Horizontal_Coordinates, Compartmentalized_Wood_Only[[k]]$Vertical_Coordinates, col = Tap_Hole_Wound_Color)
    if (Has_the_First_Tap_Been_Installed_Yet == TRUE) {
      if (k == 3) {
        break
      }
      if (k == 2) {
        Point_in_Polygon_Algorithm_1 <- as.data.frame(sapply(seq_len(nrow(Exposed_Tissue_Inside_the_Tap_Hole[[k]])), function (x) {
          sapply(seq_len((length(Tap_Hole_Wound_Polygon_Coordinates)) - 1), function (y) {
            sp::point.in.polygon(Exposed_Tissue_Inside_the_Tap_Hole[[k]][x, 1], Exposed_Tissue_Inside_the_Tap_Hole[[k]][x, 2], Tap_Hole_Wound_Polygon_Coordinates[[y]]$Horizontal_Coordinates, Tap_Hole_Wound_Polygon_Coordinates[[y]]$Vertical_Coordinates)
          })
        }))
      } else if (k > 2) {
        Point_in_Polygon_Algorithm_1 <- as.data.frame(t(sapply(seq_len(nrow(Exposed_Tissue_Inside_the_Tap_Hole[[k]])), function (x) {
          sapply(seq_len((length(Tap_Hole_Wound_Polygon_Coordinates)) - 1), function (y) {
            sp::point.in.polygon(Exposed_Tissue_Inside_the_Tap_Hole[[k]][x, 1], Exposed_Tissue_Inside_the_Tap_Hole[[k]][x, 2], Tap_Hole_Wound_Polygon_Coordinates[[y]]$Horizontal_Coordinates, Tap_Hole_Wound_Polygon_Coordinates[[y]]$Vertical_Coordinates)
          })
        })))
      }
      colnames(Point_in_Polygon_Algorithm_1) <- paste("Polygon", seq_len(ncol(Point_in_Polygon_Algorithm_1)), sep = "_")
      Overlapping_Polygon_Number_1 <- setNames(which(sapply(Point_in_Polygon_Algorithm_1, function (x) {
        any(x == 1)
      })), NULL)
      Point_in_Polygon_Algorithm_2 <- lapply(seq_len((length(Tap_Hole_Wound_Polygon_Coordinates)) - 1), function (x) {
        sapply(seq_len(nrow(Tap_Hole_Wound_Polygon_Coordinates[[x]])), function (y) {
          sp::point.in.polygon(Tap_Hole_Wound_Polygon_Coordinates[[x]][y, 1], Tap_Hole_Wound_Polygon_Coordinates[[x]][y, 2], Exposed_Tissue_Inside_the_Tap_Hole[[k]]$Horizontal_Coordinates, Exposed_Tissue_Inside_the_Tap_Hole[[k]]$Vertical_Coordinates)
        })
      })
      names(Point_in_Polygon_Algorithm_2) <- paste("Polygon", seq_len(length(Point_in_Polygon_Algorithm_2)), sep = "_")
      Overlapping_Polygon_Number_2 <- which(sapply(Point_in_Polygon_Algorithm_2, function (x) {
        any(x == 1)
      }))
      Overlapping_Polygon_Numbers <- unique(c(Overlapping_Polygon_Number_1, Overlapping_Polygon_Number_2))
      if (length(Overlapping_Polygon_Numbers) == 0) {
        Data_Frame[i, ]$Area_of_Overlap_Between_This_Tap_Hole_and_Other_Nonconductive_Wood_Columns <- 0
        Data_Frame[i, ]$Proportion_of_Exposed_Tissue_in_the_Tap_Hole_That_Overlaps_With_Other_Nonconductive_Wood_Columns <- 0
      } else if (length(Overlapping_Polygon_Numbers) > 0) {
        Data_Frame[i, ]$Area_of_Overlap_Between_This_Tap_Hole_and_Other_Nonconductive_Wood_Columns <- Calculating_Areas_of_Overlap_Function(Polygon_of_Interest = Exposed_Tissue_Inside_the_Tap_Hole[[k]], Other_Polygons = Tap_Hole_Wound_Polygon_Coordinates[Overlapping_Polygon_Numbers])
        Data_Frame[i, ]$Proportion_of_Exposed_Tissue_in_the_Tap_Hole_That_Overlaps_With_Other_Nonconductive_Wood_Columns <- Data_Frame[i, ]$Area_of_Overlap_Between_This_Tap_Hole_and_Other_Nonconductive_Wood_Columns / Cross_Sectional_Area_of_the_Exposed_Tissue_Inside_Tap_Holes
      }
    } else if (Has_the_First_Tap_Been_Installed_Yet == FALSE) {
      Data_Frame[i, ]$Area_of_Overlap_Between_This_Tap_Hole_and_Other_Nonconductive_Wood_Columns <- 0
      Data_Frame[i, ]$Proportion_of_Exposed_Tissue_in_the_Tap_Hole_That_Overlaps_With_Other_Nonconductive_Wood_Columns <- 0
      Has_the_First_Tap_Been_Installed_Yet <- TRUE
    }
  }
  i <- i + 1
}
if (Where_the_Figure_is_Being_Made == "Hard Drive") {
  dev.off()
}
write.csv(Data_Frame, file = paste0(Working_Directory, "/Sustainable Tapping Model Data.csv"))

