\name{create.plan.gate}
\alias{create.plan.gate}
\title{Generate Plan template (Gate)}
\usage{
create.plan.gate(name = "gate.simulation", ctFile.gate = NULL,
  ctWater = NULL, materialDatabase.gate = "./data/MyMaterials.db",
  HUToMaterialFile = "data/ct-HU2mat.txt", origin.gate = c(0, 0, 0),
  waterIonisationPotential.gate = 78,
  computingValues.gate = "DoseToMaterial[Gy]", saveEveryNSeconds.gate = 200,
  size.gate = "auto", position.gate = c(0, 0, 0), voxelSize.gate = c(2, 2,
  2), beams = NULL, beams.gate = "./data/beams.gate",
  particleType.gate = "proton",
  sourceDescriptionFile.gate = "./data/ProtonSimple.txt",
  totalNumberOfPrimaries = 0)
}
\description{
Generate a default (empty) template of the Plan, to be used
with Gate.
}
\seealso{
Other PlanGate: \code{\link{create.gate.structure}};
\code{\link{create.waterbox.gate}};
\code{\link{file.variable.gate}}; \code{\link{set.ct.gate}}
}

