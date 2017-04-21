library(tsensembler)

load("Data/tseries-ecml2017.rdata")
mcreps <- 10
sztr <- .5
szts <- .3
K <- 7 # 15

BL <- c('MARS', 'GLM','FFNN', 'RandomForest','PPR',
        'SVM', 'Cubist', 'GBM', 'GP')

BL.pars <- list(ffnn = list(size = c(10, 30), decay = c(0, 0.05), maxit = 750),
                mars = list(nk   = c(5, 3),   degree= c(3, 4)),
                rf = list(num.trees = c(1000, 1500), mtry = c(2, 3)),
                glm = list(alpha = c(1, 0, 0.4, 0.6)),
                gbm = list(shrinkage = c(0.005, 0.01), n.trees = c(750, 1000)),
                ppr = list(nterms = c(3, 4), sm.method = c("supsmu", "gcvspline")),
                cubist = list(committees = c(50, 100), neighbors = c(1, 3)),
                gp = list(kernel = c("rbfdot", "polydot","vanilladot"), tol = c(0.01, 0.001)),
                svm = list(kernel = c("rbfdot","polydot","vanilladot"), C = c(1, 5)))

workflow = c(Workflow(wf = 'ADE',
                      learner = BL,
                      learner.pars = BL.pars),
             Workflow(wf = 'Stacking',
                      learner = BL,
                      learner.pars = BL.pars),
             Workflow(wf = 'TSE',
                      learner = BL,
                      learner.pars = BL.pars,
                      varying.embed = FALSE,
                      varying.trainwindow = FALSE,
                      committee.ratio = NULL,
                      aggregationFUN = "static-s",
                      verbose = FALSE,
                      modelOutput = FALSE),
             Workflow(wf = 'tse.arima'),
             Workflow(wf = 'ADE_Arb',
                      learner = BL,
                      learner.pars = BL.pars),
             Workflow(wf = 'ADE_all_models',
                      learner = BL,
                      learner.pars = BL.pars),
             Workflow(wf = 'ADE_meta_runtime',
                      learner = BL,
                      learner.pars = BL.pars),
             Workflow(wf = 'Arbitrating',
                      learner = BL,
                      learner.pars = BL.pars),
             Workflow(wf = 'ADE_linear_committee',
                      learner = BL,
                      learner.pars = BL.pars),
             Workflow(wf = 'AEC',
                      learner = BL,
                      learner.pars = BL.pars),
             Workflow(wf = 'S_W50',
                      learner = BL,
                      learner.pars = BL.pars),
             Workflow(wf = 'ERP_Timmermann',
                      learner = BL,
                      learner.pars = BL.pars),
             Workflow(wf = 's_w',
                      learner = BL,
                      learner.pars = BL.pars),
             Workflow(wf = 'TSE',
                      learner = BL,
                      learner.pars = BL.pars,
                      varying.embed = FALSE,
                      varying.trainwindow = FALSE,
                      committee.ratio = NULL,
                      aggregationFUN = "static-s",
                      verbose = FALSE,
                      modelOutput = FALSE))

tse <- tsesearch(timeseries = tseries_ecml,
                 workflow = workflow,
                 embedding.dimension = K,
                 modelOutput = FALSE,
                 nReps = mcreps,
                 szTrain = sztr,
                 szTest = szts)


res <- Multi.ensembler(tse, save.each = TRUE, save.signature = "ade_experiments")
#saved

res <- ComparisonResults(sapply(res, c))
pc <- pairedComparisons(res, "ADE")

CDdiagram.BD(pc)
