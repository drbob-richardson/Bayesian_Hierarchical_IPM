# Individual-based biological generator.

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

default_fish_vital_rates <- function() {
  list(
    survival = function(size, age, sex, quality, environment, density) {
      inv_logit(-2.2 + 0.032 * size + 0.25 * quality + 0.15 * environment)
    },
    growth_mean = function(size, age, sex, quality, environment, density) {
      pmax(size, size + 18 - 0.12 * size + 1.5 * quality + environment)
    },
    growth_sd = function(size, age, sex, quality, environment, density) {
      rep(4.5, length(size))
    },
    realized_recruitment = function(
      size, age, sex, quality, environment, density
    ) {
      is_female <- sex == "F"
      is_female * exp(-7.0 + 0.075 * size + 0.20 * quality)
    },
    recruit_size = function(n, parent_size, environment) {
      rnorm(n, mean = 43 + 0.02 * parent_size + environment, sd = 3.5)
    }
  )
}

fish_vital_rates_from_inverse_ipm <- function(parameter) {
  required <- c(
    "survival_at_80",
    "survival_slope_20",
    "growth_increment_80",
    "growth_slope_20",
    "log_growth_sd",
    "log_recruitment_at_80",
    "recruitment_slope_20",
    "recruit_mean_80",
    "recruit_mean_slope_20",
    "log_recruit_sd"
  )
  stopifnot(all(required %in% names(parameter)))
  z <- function(size) (size - 80) / 20

  list(
    survival = function(size, age, sex, quality, environment, density) {
      inv_logit(
        parameter["survival_at_80"] +
          parameter["survival_slope_20"] * z(size)
      )
    },
    growth_mean = function(size, age, sex, quality, environment, density) {
      pmax(
        size,
        size + parameter["growth_increment_80"] +
          parameter["growth_slope_20"] * z(size)
      )
    },
    growth_sd = function(size, age, sex, quality, environment, density) {
      rep(exp(parameter["log_growth_sd"]), length(size))
    },
    realized_recruitment = function(
      size, age, sex, quality, environment, density
    ) {
      # Double female output so the population-average rate equals parameter.
      (sex == "F") * 2 * exp(
        parameter["log_recruitment_at_80"] +
          parameter["recruitment_slope_20"] * z(size)
      )
    },
    recruit_size = function(n, parent_size, environment) {
      rnorm(
        n,
        mean = parameter["recruit_mean_80"] +
          parameter["recruit_mean_slope_20"] * z(parent_size),
        sd = exp(parameter["log_recruit_sd"])
      )
    }
  )
}

simulate_initial_population <- function(
  n = 500,
  min_size = 35,
  max_size = 145,
  latent_quality_sd = 0
) {
  component <- sample.int(3, n, replace = TRUE, prob = c(0.62, 0.28, 0.10))
  size <- rnorm(
    n,
    mean = c(50, 78, 108)[component],
    sd = c(7, 10, 12)[component]
  )
  size <- pmin(max_size, pmax(min_size, size))

  data.frame(
    id = seq_len(n),
    birth_year = NA_integer_,
    age = sample.int(5, n, replace = TRUE) - 1L,
    sex = sample(c("F", "M"), n, replace = TRUE),
    quality = rnorm(n, sd = latent_quality_sd),
    size = size,
    stringsAsFactors = FALSE
  )
}

simulate_population <- function(
  years = 6,
  initial_population = simulate_initial_population(),
  vital_rates = default_fish_vital_rates(),
  environment = rep(0, years - 1L),
  min_size = 25,
  max_size = 160
) {
  stopifnot(years >= 2L, length(environment) == years - 1L)

  population <- initial_population
  population$birth_year[is.na(population$birth_year)] <- -population$age[
    is.na(population$birth_year)
  ]

  next_id <- max(population$id) + 1L
  census <- vector("list", years)
  transitions <- vector("list", years - 1L)
  recruits <- vector("list", years - 1L)

  for (year_index in seq_len(years)) {
    population$year <- year_index - 1L
    census[[year_index]] <- population[
      , c("year", "id", "birth_year", "age", "sex", "quality", "size")
    ]

    if (year_index == years) {
      break
    }

    env <- environment[year_index]
    density <- nrow(population)

    survival_probability <- vital_rates$survival(
      population$size,
      population$age,
      population$sex,
      population$quality,
      env,
      density
    )
    survived <- rbinom(nrow(population), 1L, survival_probability) == 1L

    growth_mean <- vital_rates$growth_mean(
      population$size,
      population$age,
      population$sex,
      population$quality,
      env,
      density
    )
    growth_sd <- vital_rates$growth_sd(
      population$size,
      population$age,
      population$sex,
      population$quality,
      env,
      density
    )
    next_size <- rep(NA_real_, nrow(population))
    next_size[survived] <- rnorm(
      sum(survived),
      mean = growth_mean[survived],
      sd = growth_sd[survived]
    )
    next_size[survived] <- pmin(
      max_size,
      pmax(min_size, next_size[survived])
    )

    transitions[[year_index]] <- data.frame(
      year = year_index - 1L,
      id = population$id,
      from_size = population$size,
      survival_probability = survival_probability,
      survived = survived,
      to_size = next_size
    )

    recruitment_mean <- vital_rates$realized_recruitment(
      population$size,
      population$age,
      population$sex,
      population$quality,
      env,
      density
    )
    offspring_count <- rpois(nrow(population), recruitment_mean)
    parent_row <- rep(seq_len(nrow(population)), offspring_count)
    n_recruits <- length(parent_row)

    if (n_recruits > 0L) {
      quality_sd <- stats::sd(population$quality)
      if (!is.finite(quality_sd)) {
        quality_sd <- 0
      }
      recruit_size <- vital_rates$recruit_size(
        n_recruits,
        population$size[parent_row],
        env
      )
      recruit_size <- pmin(max_size, pmax(min_size, recruit_size))
      new_recruits <- data.frame(
        id = seq.int(next_id, length.out = n_recruits),
        birth_year = year_index,
        age = 0L,
        sex = sample(c("F", "M"), n_recruits, replace = TRUE),
        quality = rnorm(n_recruits, sd = quality_sd),
        size = recruit_size,
        stringsAsFactors = FALSE
      )
      recruits[[year_index]] <- data.frame(
        year = year_index,
        id = new_recruits$id,
        parent_id = population$id[parent_row],
        parent_size = population$size[parent_row],
        size = new_recruits$size
      )
      next_id <- next_id + n_recruits
    } else {
      new_recruits <- population[FALSE, c(
        "id", "birth_year", "age", "sex", "quality", "size"
      )]
      recruits[[year_index]] <- data.frame(
        year = integer(),
        id = integer(),
        parent_id = integer(),
        parent_size = numeric(),
        size = numeric()
      )
    }

    survivors <- population[survived, c(
      "id", "birth_year", "age", "sex", "quality", "size"
    )]
    survivors$age <- survivors$age + 1L
    survivors$size <- next_size[survived]
    population <- rbind(survivors, new_recruits)
  }

  list(
    census = do.call(rbind, census),
    transitions = do.call(rbind, transitions),
    recruits = do.call(rbind, recruits),
    environment = environment
  )
}
