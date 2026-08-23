#' Fit an XGBoost diagnostic classifier and rank features
#'
#' Select a learning rate using cross-validated binary log loss, fit a final
#' XGBoost classifier, and create a feature-gain importance plot.
#'
#' @param data A data frame containing predictor columns and one categorical
#'   response column.
#' @param feature_colnames Character vector naming predictor features used by
#'   XGBoost.
#' @param response_colnames Name of the categorical response column. Although
#'   the argument name is plural, the current implementation expects one name.
#' @param n_folds Number of cross-validation folds passed to
#'   [xgboost::xgb.cv()].
#' @param seed Integer random seed set before cross-validation and model
#'   fitting.
#'
#' @return A named list containing `result`, an XGBoost feature-importance table
#'   sorted by decreasing gain; `plot`, a horizontal gain bar plot; `best_eta`,
#'   the selected learning rate; and `best_cv_logloss`, the smallest recorded
#'   cross-validated test log loss for that learning rate.
#'
#' @details
#' Response levels are ordered with `group_new_levels()` and converted to
#' binary labels: the first level becomes 0 and every other level becomes 1.
#' The function is intended for binary classification. Predictor data are
#' converted to a double matrix, and missing or other non-finite matrix values
#' are replaced by zero.
#'
#' Four learning rates (`0.001`, `0.01`, `0.1`, and `0.3`) are evaluated with
#' 30-round cross-validation. For each rate, the minimum test log loss across
#' rounds is used as its score. Other tuning parameters are fixed: maximum tree
#' depth 6, minimum child weight 1, subsample 0.5, column subsample 1, and gamma
#' 0. The final binary logistic model uses the selected learning rate and 200
#' boosting rounds.
#'
#' Feature importance is calculated with [xgboost::xgb.importance()] and sorted
#' by `Gain`. Periods in the feature names prepared for the importance table are
#' displayed as hyphens. Predictor columns should therefore be numeric and
#' their model-matrix representation should correspond to the matrix used to
#' train the model.
#'
#' @examples
#' \dontrun{
#' diagnosis_data <- data.frame(
#'   GeneA = rnorm(80),
#'   GeneB = rnorm(80),
#'   GeneC = rnorm(80),
#'   Group = factor(rep(c("Control", "Disease"), each = 40))
#' )
#'
#' xgb_result <- xgboost_diagnosis(
#'   data = diagnosis_data,
#'   feature_colnames = c("GeneA", "GeneB", "GeneC"),
#'   response_colnames = "Group",
#'   n_folds = 5,
#'   seed = 123
#' )
#'
#' xgb_result$best_eta
#' xgb_result$best_cv_logloss
#' xgb_result$result
#' xgb_result$plot
#' }
#'
#' @seealso [xgboost::xgb.cv()], [xgboost::xgb.train()],
#'   [xgboost::xgb.importance()]
#' @export
xgboost_diagnosis <- function(data = df_train, feature_colnames, response_colnames,
                              n_folds = 5, seed = 1){

  suppressPackageStartupMessages({
    library(dplyr)
    library(xgboost)
    library(Matrix)
    library(ggplot2)
  })

  # data process（沿用你的逻辑）
  data[[response_colnames]] <- factor(
    data[[response_colnames]],
    levels = group_new_levels(vector = as.character(data[[response_colnames]]))
  )
  data[[response_colnames]] <- ifelse(
    data[[response_colnames]] == group_new_levels(vector = as.character(data[[response_colnames]]))[1],
    0, 1
  )
  data[[response_colnames]] <- factor(data[[response_colnames]], levels = c(0,1))

  set.seed(seed)

  # grid（沿用你的逻辑：只调 eta，nrounds 固定 30）
  etas <- c(0.001, 0.01, 0.1, 0.3)
  nrounds_tune <- 30

  # ==== 用 xgb.cv 选择最佳 eta（替代 caret::train）====
  X_cv <- as.matrix(data[, feature_colnames, drop = FALSE])
  storage.mode(X_cv) <- "double"
  X_cv[!is.finite(X_cv)] <- 0

  y_cv <- as.numeric(data[, response_colnames]) - 1
  stopifnot(all(y_cv %in% c(0,1)))

  dcv <- xgb.DMatrix(data = X_cv, label = y_cv)

  best_eta <- etas[1]
  best_score <- Inf

  for (e in etas) {
    cv <- xgb.cv(
      params = list(
        objective = "binary:logistic",
        eval_metric = "logloss",
        eta = e,
        gamma = 0,
        max_depth = 6,
        min_child_weight = 1,
        subsample = 0.5,
        colsample_bytree = 1
      ),
      data = dcv,
      nrounds = nrounds_tune,
      nfold = n_folds,
      verbose = 0
    )

    log <- cv$evaluation_log
    col <- intersect(c("test_logloss_mean", "test_logloss"), colnames(log))[1]
    score <- min(log[[col]], na.rm = TRUE)

    if (is.finite(score) && score < best_score) {
      best_score <- score
      best_eta <- e
    }
  }

  # ==== 训练最终模型（沿用你原来 200 轮的思路）====
  traindata1 <- data.matrix(data[, feature_colnames, drop = FALSE])
  storage.mode(traindata1) <- "double"
  traindata1[!is.finite(traindata1)] <- 0
  traindata2 <- Matrix(traindata1, sparse = TRUE)

  train_y <- as.numeric(data[, response_colnames]) - 1
  dtrain <- xgb.DMatrix(data = traindata2, label = train_y)

  params <- list(
    objective = "binary:logistic",
    eta = best_eta,
    colsample_bytree = 1,
    min_child_weight = 1,
    gamma = 0,
    subsample = 0.5,
    max_depth = 6
  )

  res.xgb <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 200,
    verbose = 0
  )

  # importance：保留你“把 . 换成 - ”的呈现逻辑
  ddf_reain <- data
  colnames(ddf_reain) <- make.names(colnames(ddf_reain))
  train_matrix <- sparse.model.matrix(as.formula(paste(response_colnames, "~ . -1")), data = ddf_reain)
  feature_names_for_importance <- gsub("[.]", "-", train_matrix@Dimnames[[2]])

  xgb_importance_df <- xgb.importance(
    feature_names = feature_names_for_importance,
    model = res.xgb
  )

  xgb_plot_data <- xgb_importance_df %>% dplyr::arrange(desc(Gain))

  ggplot(xgb_plot_data, aes(x = reorder(Feature, Gain), y = Gain)) +
    geom_bar(aes(fill = Gain), stat = "identity") +
    scale_fill_viridis_c() +
    coord_flip() +
    labs(title = "XGBoost Feature Importance", x = "Gene", y = "Importance") +
    theme_classic(base_size = 16) +
    theme(axis.title = element_text(face = "bold"),
          legend.position = "none",
          plot.margin = ggplot2::margin(t = 10, r = 20, b = 10, l = 10, unit = "pt")
         ) -> p

  return(list(result = xgb_plot_data, plot = p, best_eta = best_eta, best_cv_logloss = best_score))
}
