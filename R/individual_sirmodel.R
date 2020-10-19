#' @importFrom stats rbinom
bernoulli <- function(size, p) sample.int(size, rbinom(1, size, p))


#' Stochastic SIR model, individual
#'
#' Takes timestep steps over time for population N
#'
#' @param end_time end time for data
#' @param pars parameter list
#'
#' @return dataframe
#'
#' @export
#'
#' @examples
#' individual_sirmodel(100, list())
individual_sirmodel <- function(end_time, pars = NULL, vars = NULL, events = NULL) {


  infection_process <- function(api) {

    pars <- api$get_parameters()
    FOI <- pars$beta * length(api$get_state(human, I)) / pars$N
    susceptible <- api$get_state(human, S)
    infected <- bernoulli(length(susceptible), 1.0 - exp(-1 * FOI * pars$dt))
    api$queue_state_update(human, I, susceptible[infected])

  }

  recovery_process <- function(api) {

    pars <- api$get_parameters()
    infected <- api$get_state(human, I)
    recovered <- bernoulli(length(infected), 1.0 - exp(-1 * pars$nu * pars$dt))
    api$queue_state_update(human, R, infected[recovered])

  }

  death_process <- function(api) {

    pars <- api$get_parameters()
    if (pars$includebirth) {
      died <- bernoulli(pars$N, 1.0 - exp(-1 * pars$mu * pars$dt))
      api$queue_state_update(human, S, died)
    }
  }

  processes <- list(infection_process, recovery_process, death_process)

  render_state_sizes <- function(api) {
    api$render('susceptable_counts', length(api$get_state(human, S)))
    api$render('infected_counts', length(api$get_state(human, I)))
    api$render('recovered_counts', length(api$get_state(human, R)))
  }


  # initialist pars
  pars <- get_parameters(pars)


  # Set no of infections, NI, no of susceptible, pops,
  # no of recovered at the start and no of time points
  NI <- pars$I0
  population <- pars$N
  pops <- population - NI
  timestep <- pars$num/pars$dt
  S <- State$new('S', pops)
  I <- State$new('I', NI)
  R <- State$new('R', 0)

  # Variables
  human <- Individual$new(name = 'human',
                          states = list(S, I, R))

  output <- simulate(individuals = human,
                     processes = list(
                       infection_process,
                       recovery_process,
                       death_process,
                       render_state_sizes
                     ),
                     end_timestep  = timestep,
                     parameters = pars)

  df <-   data.frame(S = output$susceptable_counts,
                     I = output$infected_counts,
                     R = output$recovered_counts,
                     time = output$time,
                     type = "Individual",
                     legend = "Individual",
                     stringsAsFactors = FALSE)

}
