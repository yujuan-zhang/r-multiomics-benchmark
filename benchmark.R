# Two-block classification benchmark. All learned preprocessing is training-only.
prepare_demo <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  e <- new.env(); data("nutrimouse", package = "mixOmics", envir = e)
  d <- e$nutrimouse; id <- sprintf("mouse_%02d", seq_len(nrow(d$gene)))
  write.csv(data.frame(sample_id = id, outcome = d$diet, stratum = d$genotype),
            file.path(path, "samples.csv"), row.names = FALSE)
  write.csv(data.frame(sample_id = id, d$gene, check.names = FALSE),
            file.path(path, "block_a.csv"), row.names = FALSE)
  write.csv(data.frame(sample_id = id, log1p(d$lipid), check.names = FALSE),
            file.path(path, "block_b.csv"), row.names = FALSE)
  writeLines(c("Nutrimouse, distributed with mixOmics; export version:",
    as.character(utils::packageVersion("mixOmics")),
    "https://bioconductor.org/packages/mixOmics/",
    "40 mice; 120 preselected gene-expression measurements; 21 hepatic fatty acids.",
    "Outcome: five diets. Fold-balancing stratum: genotype (wt / ppar).",
    "Block A: supplied expression values from nylon macroarrays (not RNA-seq).",
    "Block B: log1p of supplied fatty-acid percentages; zeros retained as zero.",
    "This is not a log-ratio analysis; compositional dependence remains.",
    "IDs mouse_01..mouse_40 preserve the original package row order.",
    "Martin et al. (2007), Novel aspects of PPARalpha-mediated regulation of lipid",
    "and xenobiotic metabolism revealed through a multigenomic study.",
    "Demo reproduces a historical, preprocessed dataset; preprocessing upstream",
    "of the published matrices cannot be repeated within cross-validation."),
    file.path(path, "PROVENANCE.txt"))
}

load_input <- function(path) {
  s <- read.csv(file.path(path, "samples.csv"), stringsAsFactors = FALSE)
  if (!all(c("sample_id", "outcome", "stratum") %in% names(s)) ||
      anyNA(s) || anyDuplicated(s$sample_id) || any(!nzchar(s$sample_id)))
    stop("samples.csv requires unique nonempty sample_id, outcome and stratum.")
  x <- lapply(c(a = "block_a.csv", b = "block_b.csv"), function(f) {
    z <- read.csv(file.path(path, f), check.names = FALSE)
    if (!identical(names(z)[1], "sample_id") || anyDuplicated(z[[1]]) ||
        !setequal(z[[1]], s$sample_id) || anyDuplicated(names(z)[-1]))
      stop("Each block must contain unique features and exactly the metadata sample IDs.")
    z <- z[match(s$sample_id, z[[1]]), -1, drop = FALSE]
    if (ncol(z) < 2 || !all(vapply(z, is.numeric, logical(1))))
      stop("Each block needs at least two numeric features.")
    z <- as.matrix(z); rownames(z) <- s$sample_id
    if (any(!is.finite(z))) stop("Missing/non-finite values are unsupported; supply complete paired matrices.")
    z
  })
  s$outcome <- factor(s$outcome)
  if (nlevels(s$outcome) < 2) stop("At least two outcome classes are required.")
  list(x = x, s = s)
}
folds <- function(y, stratum, k, seed) {
  set.seed(seed); group <- interaction(y, stratum, drop = TRUE)
  if (any(table(group) < k)) stop("Each outcome x stratum cell needs at least k samples.")
  out <- integer(length(y))
  for (g in levels(group)) {
    ii <- which(group == g); out[ii] <- sample(rep(seq_len(k), length.out = length(ii)))
  }
  out
}
preprocess <- function(train, test) {
  mu <- colMeans(train); sd <- apply(train, 2, stats::sd)
  keep <- is.finite(sd) & sd > 1e-8
  if (sum(keep) < 2) stop("Fewer than two variable training features in a block.")
  transform <- function(z) sweep(sweep(z[, keep, drop = FALSE], 2, mu[keep]), 2, sd[keep], "/")
  list(train = transform(train), test = transform(test))
}
partition <- function(x, tr, te) {
  p <- lapply(x, function(z) preprocess(z[tr, , drop = FALSE], z[te, , drop = FALSE]))
  list(train = lapply(p, `[[`, "train"), test = lapply(p, `[[`, "test"))
}
score <- function(y, pred) {
  if (anyNA(pred) || !all(pred %in% levels(y))) stop("Invalid predicted classes.")
  mean(vapply(levels(y), function(k) mean(pred[y == k] == k), numeric(1)))
}
joined <- function(x, method) {
  if (method == "A only") return(x$a)
  if (method == "B only") return(x$b)
  z <- do.call(cbind, lapply(names(x), function(b) {
    m <- x[[b]]; colnames(m) <- paste0(b, "::", colnames(m)); m
  })); z
}
fit_enet <- function(x, y, alpha, lambda) {
  withCallingHandlers(glmnet::glmnet(x, y, family = if (nlevels(y) == 2) "binomial" else "multinomial",
    alpha = alpha, lambda = sort(unique(c(10, lambda)), decreasing = TRUE),
    standardize = FALSE, control = list(thresh = 1e-7, maxit = 100000)),
    warning = function(w) {
      # Announced once by run_benchmark and retained in the report.
      if (grepl("class has fewer than 8", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    })
}
fit_diablo <- function(x, y, grid) {
  design <- matrix(grid$design, 2, 2, dimnames = list(names(x), names(x))); diag(design) <- 0
  keep <- lapply(x, function(z) rep(min(grid$budget, ncol(z)), grid$ncomp))
  suppressMessages(mixOmics::block.splsda(x, y, ncomp = grid$ncomp,
    keepX = keep, design = design, scale = FALSE, max.iter = 1000))
}
predict_diablo <- function(model, x, ncomp) {
  p <- withCallingHandlers(predict(model, newdata = x), warning = function(w) {
    if (startsWith(conditionMessage(w), "Some blocks are missing")) invokeRestart("muffleWarning")
  })
  votes <- if (length(x) == 1) p$class else p$WeightedVote
  as.character(votes$centroids.dist[, ncomp])
}

run_benchmark <- function(input, output, cfg) {
  for (pkg in c("mixOmics", "glmnet", "ggplot2", "patchwork", "ragg"))
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install dependency: ", pkg)
  if (utils::packageVersion("glmnet") < "5.1") stop("glmnet >= 5.1 is required for the control API.")
  options(rgl.useNULL = TRUE)
  d <- load_input(input); x <- d$x; s <- d$s; y <- s$outcome
  message("Small classes can make glmnet estimates unstable; repeated <8-observation class warnings are consolidated here. See README limitations.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  eg <- expand.grid(alpha = cfg$alpha, lambda = cfg$lambda)
  dg <- expand.grid(budget = cfg$budgets, ncomp = cfg$components, design = cfg$designs)
  predictions <- selections <- tuning <- splits <- list()
  add_pred <- function(pred, method, rep, fold, te, count = NA_integer_) {
    predictions[[length(predictions) + 1L]] <<- data.frame(repeat_id = rep, fold = fold,
      sample_id = s$sample_id[te], truth = as.character(y[te]), method = method,
      prediction = pred, selected_features = count)
  }
  for (rr in seq_len(cfg$repeats)) {
    outer <- folds(y, s$stratum, cfg$outer_folds, cfg$seed + rr)
    for (ff in seq_len(cfg$outer_folds)) {
      message(sprintf("Repeat %d/%d | outer fold %d/%d", rr, cfg$repeats, ff, cfg$outer_folds))
      tr <- which(outer != ff); te <- which(outer == ff)
      inner <- folds(y[tr], s$stratum[tr], cfg$inner_folds, cfg$seed + 100L * rr + ff)
      splits[[length(splits) + 1]] <- data.frame(repeat_id = rr, outer_fold = ff,
        sample_id = s$sample_id, role = ifelse(outer == ff, "test", "train"),
        inner_fold = replace(rep(NA_integer_, nrow(s)), tr, inner))
      pp <- lapply(seq_len(cfg$inner_folds), function(f) {
        a <- tr[inner != f]; b <- tr[inner == f]
        c(partition(x, a, b), list(y = y[a], valid_y = y[b], idx = which(inner == f)))
      })
      full <- partition(x, tr, te)
      best_scores <- numeric(3); names(best_scores) <- c("A only", "B only", "Concatenated")
      final_predictions <- list()
      for (method in names(best_scores)) {
        scores <- vapply(seq_len(nrow(eg)), function(g) {
          pred <- character(length(tr))
          for (p in pp) {
            m <- fit_enet(joined(p$train, method), p$y, eg$alpha[g], eg$lambda[g])
            pred[p$idx] <- as.character(predict(m, joined(p$test, method), s = eg$lambda[g], type = "class"))
          }
          score(y[tr], pred)
        }, numeric(1))
        # Ties prefer stronger regularization, then smaller alpha (predeclared grid order).
        win <- which(scores == max(scores)); win <- win[order(-eg$lambda[win], eg$alpha[win])][1]
        best_scores[method] <- scores[win]
        tuning[[length(tuning) + 1]] <- data.frame(repeat_id = rr, fold = ff, method = method,
          parameter = sprintf("alpha=%s;lambda=%s", eg$alpha, eg$lambda),
          inner_balanced_accuracy = scores, chosen = seq_len(nrow(eg)) == win)
        model <- fit_enet(joined(full$train, method), y[tr], eg$alpha[win], eg$lambda[win])
        pr <- as.character(predict(model, joined(full$test, method), s = eg$lambda[win], type = "class"))
        final_predictions[[method]] <- pr
        add_pred(pr, method, rr, ff, te)
      }
      chosen_single <- names(which.max(best_scores[c("A only", "B only")]))
      tuning[[length(tuning) + 1]] <- data.frame(repeat_id = rr, fold = ff,
        method = "Inner-selected single", parameter = c("A only", "B only"),
        inner_balanced_accuracy = best_scores[c("A only", "B only")],
        chosen = c("A only", "B only") == chosen_single)
      add_pred(final_predictions[[chosen_single]], "Inner-selected single", rr, ff, te)
      scores <- vapply(seq_len(nrow(dg)), function(g) {
        pred <- character(length(tr))
        for (p in pp) {
          m <- fit_diablo(p$train, p$y, dg[g, ])
          pred[p$idx] <- predict_diablo(m, p$test, dg$ncomp[g])
        }
        score(y[tr], pred)
      }, numeric(1))
      # Evaluate each predeclared budget separately and a budget-tuned primary model.
      for (budget in c(NA, cfg$budgets)) {
        eligible <- if (is.na(budget)) seq_len(nrow(dg)) else which(dg$budget == budget)
        win <- eligible[which.max(scores[eligible])]
        method <- if (is.na(budget)) "DIABLO" else paste0("Panel ", budget, "/block/component")
        tuning[[length(tuning) + 1]] <- data.frame(repeat_id = rr, fold = ff, method = method,
          parameter = sprintf("keep=%s;ncomp=%s;design=%s", dg$budget, dg$ncomp, dg$design),
          inner_balanced_accuracy = scores, chosen = seq_len(nrow(dg)) == win)
        model <- fit_diablo(full$train, y[tr], dg[win, ])
        features <- do.call(rbind, lapply(names(x), function(b) {
          load <- model$loadings[[b]]
          selected <- rownames(load)[rowSums(abs(load) > 1e-10) > 0]
          data.frame(block = b, feature = selected)
        }))
        selections[[length(selections) + 1]] <- cbind(repeat_id = rr, fold = ff, method, features)
        add_pred(predict_diablo(model, full$test, dg$ncomp[win]), method, rr, ff, te, nrow(features))
        if (is.na(budget)) for (b in names(x))
          add_pred(predict_diablo(model, full$test[b], dg$ncomp[win]),
            paste0("DIABLO: ", toupper(b), " available"), rr, ff, te, nrow(features))
      }
    }
  }
  pred <- do.call(rbind, predictions); sel <- do.call(rbind, selections)
  metrics <- do.call(rbind, lapply(split(pred, interaction(pred$repeat_id, pred$method, drop = TRUE)), function(z)
    data.frame(repeat_id = z$repeat_id[1], method = z$method[1],
      balanced_accuracy = score(factor(z$truth, levels = levels(y)), z$prediction),
      mean_selected_features = if (all(is.na(z$selected_features))) NA_real_ else mean(z$selected_features))))
  write.csv(metrics, file.path(output, "metrics.csv"), row.names = FALSE)
  for (item in c("predictions", "selected_features", "tuning", "splits")) {
    z <- switch(item, predictions = pred, selected_features = sel,
      tuning = do.call(rbind, tuning), splits = do.call(rbind, splits))
    con <- gzfile(file.path(output, paste0(item, ".csv.gz")), "wt")
    write.csv(z, con, row.names = FALSE); close(con)
  }
  overview(metrics, sel, cfg, file.path(output, "overview.png"))
  means <- aggregate(balanced_accuracy ~ method, metrics, mean)
  gain <- means$balanced_accuracy[means$method == "DIABLO"] -
    means$balanced_accuracy[means$method == "Inner-selected single"]
  writeLines(c("# Benchmark results", "", sprintf("%d samples; %d classes; %d x %d-fold nested cross-validation.",
    nrow(s), nlevels(y), cfg$repeats, cfg$outer_folds), "",
    "| Method | Mean balanced accuracy |", "|---|---:|",
    sprintf("| %s | %.3f |", means$method, means$balanced_accuracy), "",
    sprintf("DIABLO minus inner-selected single-omics baseline: %+.1f percentage points.", 100 * gain),
    "", "Repeated estimates share samples and training sets; their spread is descriptive, not a confidence interval.",
    "Panel budgets apply per block per component. Actual panel size is the union across components.",
    "Missing-block scenarios reuse the complete-data-trained DIABLO model without retraining.",
    "Selection frequency measures model reuse across overlapping training sets, not biomarker validation.",
    "Small training classes (<8 observations in this demo) can make glmnet estimates unstable.",
    "These results describe this dataset and search grid; no method is assumed to win."), file.path(output, "summary.md"))
  capture.output(list(configuration = cfg, input_md5 = tools::md5sum(file.path(input,
    c("samples.csv", "block_a.csv", "block_b.csv"))), session = sessionInfo()),
    file = file.path(output, "run_info.txt"))
  message(sprintf("Complete. DIABLO gain over inner-selected single: %+.1f percentage points.", 100 * gain))
  invisible(metrics)
}

overview <- function(metrics, selections, cfg, path) {
  library(ggplot2)
  theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"), plot.title.position = "plot"))
  primary <- subset(metrics, method %in% c("A only", "B only", "Concatenated", "DIABLO", "Inner-selected single"))
  primary$method <- factor(primary$method, levels = c("A only", "B only", "Concatenated", "Inner-selected single", "DIABLO"))
  a <- ggplot(primary, aes(method, balanced_accuracy, group = repeat_id)) +
    geom_line(color = "#B8C7D3", linewidth = .6) + geom_point(color = "#176B87", size = 2.6) +
    coord_flip(ylim = c(0, 1)) + labs(title = "Does integration improve prediction?", x = NULL, y = "Balanced accuracy")
  missing <- subset(metrics, grepl("DIABLO", method))
  b <- ggplot(missing, aes(method, balanced_accuracy, group = repeat_id)) +
    geom_line(color = "#D6C7B7") + geom_point(color = "#C16E35", size = 2.6) +
    coord_flip(ylim = c(0, 1)) + labs(title = "What happens when an assay is unavailable?", x = NULL, y = "Balanced accuracy")
  panels <- subset(metrics, grepl("Panel", method))
  c <- ggplot(panels, aes(mean_selected_features, balanced_accuracy, color = method)) +
    geom_point(size = 3) + scale_color_manual(values = c("#176B87", "#C16E35")) +
    coord_cartesian(ylim = c(0, 1)) + labs(title = "Panel size versus performance", x = "Mean unique selected features (both blocks)",
      y = "Balanced accuracy", color = NULL) + theme(legend.position = "bottom", legend.text = element_text(size = 8))
  sel <- subset(selections, method == "DIABLO")
  freq <- aggregate(fold ~ block + feature, sel, length)
  freq$frequency <- freq$fold / (cfg$repeats * cfg$outer_folds)
  freq <- head(freq[order(-freq$frequency, freq$block, freq$feature), ], 12)
  freq$label <- factor(paste(freq$block, freq$feature, sep = ": "), levels = rev(paste(freq$block, freq$feature, sep = ": ")))
  d <- ggplot(freq, aes(frequency, label, fill = block)) + geom_col(width = .7) +
    scale_fill_manual(values = c("#176B87", "#C16E35")) + scale_x_continuous(limits = c(0, 1)) +
    labs(title = "Which features recur across training sets?", x = "Selection frequency · primary DIABLO", y = NULL) +
    theme(legend.position = "none")
  p <- patchwork::wrap_plots(a, b, c, d, ncol = 2) + patchwork::plot_annotation(
    title = "Multi-omics: added value, resilience and panel size",
    subtitle = "Nested cross-validation | paired outer splits | held-out predictions only",
    caption = "A: expression · B: lipids in the bundled demo. Repeats overlap; points are not independent replicates.",
    tag_levels = "A")
  ggsave(path, p, width = 13, height = 9, dpi = 180, device = ragg::agg_png, bg = "white")
}

check_contracts <- function() {
  y <- factor(rep(c("a", "b"), each = 8)); z <- rep(rep(c("u", "v"), each = 4), 2)
  f <- folds(y, z, 4, 1)
  stopifnot(all(table(interaction(y, z), f) == 1), identical(f, folds(y, z, 4, 1)))
  x <- matrix(seq_len(48), 16, 3)
  p <- preprocess(x[1:12, ], x[13:16, ])
  q <- preprocess(x[1:12, ], x[13:16, ] * 1000)
  stopifnot(identical(p$train, q$train), max(abs(colMeans(p$train))) < 1e-10,
            score(y, as.character(y)) == 1, score(y, rep("a", 16)) == .5)
  bad <- try(folds(y[1:3], z[1:3], 4, 1), silent = TRUE)
  stopifnot(inherits(bad, "try-error"))
  # Verify keyed alignment, and reject duplicates/non-finite inputs.
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  s <- data.frame(sample_id = paste0("s", 1:16), outcome = y, stratum = z)
  write.csv(s, file.path(tmp, "samples.csv"), row.names = FALSE)
  block <- data.frame(sample_id = s$sample_id, x)
  write.csv(block[16:1, ], file.path(tmp, "block_a.csv"), row.names = FALSE)
  write.csv(block, file.path(tmp, "block_b.csv"), row.names = FALSE)
  loaded <- load_input(tmp)
  stopifnot(identical(unname(loaded$x$a), unname(loaded$x$b)))
  block$sample_id[2] <- block$sample_id[1]
  write.csv(block, file.path(tmp, "block_b.csv"), row.names = FALSE)
  stopifnot(inherits(try(load_input(tmp), silent = TRUE), "try-error"))
  message("Checks passed: stratification, deterministic splits, training-only scaling, metrics and sample alignment.")
}
