
cleanWindRaw <- function(raw.data, col=c('Close'), ...) {
# Clean raw dataset from wind TDB api.
#
  h.digits = 10e6
  m.digits = 10e4
  raw.data$Hour = (raw.data$Time) %/% h.digits
  raw.data$Min = (raw.data$Time - h.digits * raw.data$Hour) %/% m.digits
  # Export timeseries.
  wind.xts = xts(raw.data[, col], as.POSIXct(
  					paste(raw.data$Date,
  						  raw.data$Hour,
  						  raw.data$Min),
  					format = "%Y%m%d %H %M"))
  wind.xts
}

windM1Data <- function(...) UseMethod('windM1Data')

windM1Data.default <- function(...) {
# Export a list of data, ts and params for modeling.
#
  C1 = cleanWindRaw(rawIFC1) / 10e3
  C2 = cleanWindRaw(rawIFC2) / 10e3 # Bp -> p
  HS300 = cleanWindRaw(rawHS300) / 10e3
  dC = C2 - C1 # next month - current month
  basis = C1 - HS300 # future - underlying
  merged = merge(C1, C2, 'inner')
  set = list(
  	C1.xts = C1, C2.xts = C2,
  	HS300.xts = HS300,
  	basis = basis,
  	diff = dC, merged = merged
  )
  class(set) = 'windM1Data'
  set
}

plot.windM1Data <- function(data, which=1, ...) {
  if(which==1) {
  	IFC1 = data$C1.xts
  	chartSeries(IFC1, theme = chartTheme('white'),
  				TA=c(addTA(RSI(data$C1.xts))))
  } else if(which==2) {
  	IFC2 = data$C2.xts
  	chartSeries(IFC2, theme = chartTheme('white'),
  				TA = c(addTA(RSI(IFC2)), addTA(MACD(IFC2))))
  } else if(which==3) {

  }
}

zGARCHsignal <- function(x, start, window.size,
						  n.ahead=1, spec=NULL, cluster=NULL) {
# Params
# ------
# x: (xts object) the timeseries at which we target
# start: (int) estimation starting point, i.e. the size of
# 	in.the.sample part for the first dynamic estimation.
# window.size: (int) the number of period before re-estimating model
# n.ahead: (int) the forecast horizon.
# spec: (rugarch.spec object) model specification
# cluster: parallel cluster.
#
# Params in the Scope
# -------------------
# x: the target ts.
# N: nrow(x).
# reest.point: vector of time at which model is re-estimated.
# k: length(reest.point)
# results
# forecast
  if (!is.xts(x))
	stop('\n[garch]: must pass an xts object to x') # require x be xts.
  x = na.omit(x)
  N = nrow(x) # full sample size
  if (N < start)
  	stop('\n[garch]: sample size < starting time.')
  # re-estimate time point. (as index)
  reest.point = seq(start+window.size, N, by=window.size)
  k = length(reest.point) # re-estimation timings
  # out-of-sample data size for each estimation
  out.sample.size = rep(window.size, k)

  if(reest.point[k]<N){ # add the rest of data to re-estimation schedule
  	reest.point = c(reest.point, N)
  	k = length(reest.point)
  	out.sample.size = c(out.sample.size, reest.point[k]-reest.point[k-1])
  }

  if(is.null(spec)){ # set default model if not given.
  	spec = ugarchspec(variance.model = list(
  	  model = 'sGARCH'
  	), distribution.model = 'std')
  } # as ARMA(1,1)-sGARCH(1,1), innovation ~ Student t

  if(!is.null(cluster)){ # === parallel computation ===
  	#stop('[garch]: NotImplemented!')
  	clusterEvalQ(cluster, library(rugarch))
  	clusterEvalQ(cluster, library(xts))
  	clusterExport(cluster, varlist = c('x', 'reest.point',
  				'window.size', 'spec', 'out.sample.size'), envir = environment())
  	# ------ (CL)Main Loop ------
  	results = clusterApply(cl = cluster, 1:k, fun = function(i) {
  	  model = ugarchfit(spec,
  	  				    x[1:reest.point[i],],
  	  				    out.sample = out.sample.size[i],
  	  				    solver = 'hybrid')
  	  ahead = ugarchforecast(model, n.ahead = 1, n.roll = out.sample.size[i]-1)
  	  conditional.mean = fitted(ahead)
  	  conditional.sigma = sigma(ahead)
  	  # paste value to the out.sample positions.
  	  conditional.mean.xts = xts(as.numeric(conditional.mean), tail(
  	  	time(x[1:reest.point[i], ]), out.sample.size[i]))
  	  conditional.sigma.xts = xts(as.numeric(conditional.sigma), tail(
  	  	time(x[1:reest.point[i], ]), out.sample.size[i]))
  	  y = list(
  	  	mean = conditional.mean.xts,
  	  	sigma = conditional.sigma.xts
  	  )
  	  print(paste('[RollingModel]: Finished', i))
  	  return(y)
  	}) # ------ End (CL)Main Loop ------
  } else { # === unparalleled version ===
  	# ------ Main Loop ------
  	#
  	results = lapply(1:k, # for k re-estimation timing
  	  FUN = function(i) {
  	  	model = ugarchfit(spec,
  	  					  x[1:reest.point[i],], # re-estimate for all history.
  	  					  out.sample = out.sample.size[i], # out of sample data
  	  					  solver = 'hybrid')
  	  	ahead = ugarchforecast(model, n.ahead = 1, # forecast on out.of.sample points.
  	  						   n.roll = out.sample.size[i] - 1)
  	  	conditional.mean = fitted(ahead)
  	  	conditional.sigma = sigma(ahead)
  	  	# paste value to the out.sample positions.
  	  	conditional.mean.xts = xts(as.numeric(conditional.mean), tail(
  	  		time(x[1:reest.point[i], ]), out.sample.size[i]))
  	  	conditional.sigma.xts = xts(as.numeric(conditional.sigma), tail(
  	  		time(x[1:reest.point[i], ]), out.sample.size[i]))
  	  	y = list(
  	  	  mean = conditional.mean.xts,
  	  	  sigma = conditional.sigma.xts
  	  	)
  	  	print(paste('[RollingModel]: Finished', i))
  	  	return(y)
  	  })
  	# ------ End Main Loop ------
  } # end if

  # Export modeling results
  forecast = merge(results[[1]]$mean, results[[1]]$sigma)
  for (i in 2:k) {
  	forecast = rbind(forecast, merge(results[[i]]$mean, results[[i]]$sigma))
  }
  colnames(forecast) = c('cond.mean', 'cond.std')
  zg.signal = (x - forecast$cond.mean) / forecast$cond.std
  zgarch.results = list(
  	forecast = forecast,
  	signal = zg.signal,
  	spec = spec,
  	start = start,
  	roll = window.size,
  	n.reest = k,
  	n.sample = N,
  	x = x
  )
  class(zgarch.results) = 'zgarch'
  return(zgarch.results)
}


plotSignal.zgarch = function(zgarch, y.qtl=0.999, thres=c(0.5, 2),...) {
	# plot Arb trading signals.
	# Params
	# ------
	# * zgarch: (zgarch object) estimation result of rolling model.
	# * thres: (vector of numerics) trading rule threshold, in which
	# 	- +-thres[1]: signal of sell/cover, default 0.5
	# 	- +-thres[2]: signal of long/short, default 2
	# means that Long order when series score (X-\hat{E}(X))/\hat{std}(X)
	# breaks -thres[2], i.e. X deviates thres[2] * sigma from its cond.mean
	# * y.qtl: (float) a confidence level, whose corresponding quantile is
	# 	treated as ylim in plots.
	ord.index = order(index(zgarch$forecast)) # order indices
	demean = zgarch$x - zgarch$forecast$cond.mean
	std = zgarch$forecast$cond.std
	score = demean / std
	y.max = quantile(as.vector(thres[2]*std), y.qtl)
	y.lim = c(-y.max, y.max)
	x.lim = c(0, length(ord.index))
	# plot demean + +-n*std
	plot0 = xyplot(demean~ord.index,
				   ylim=y.lim, xlim=x.lim,
				   col='black',
				   xlab='time(periodCounts)',
				   ylab='value',
				   main='Arbitrage Signals',
				   type='l')
	plot.sig.up = xyplot(thres[2]*std~ord.index,
						 ylim=y.lim, xlim=x.lim,
						 type='l')
	plot.sig.down = xyplot(-thres[2]*std~ord.index,
						   ylim=y.lim, xlim=x.lim,
						   type='l')
	# plot signals
	tmp = data.frame(ord.index, demean, score)
	colnames(tmp) = c('ord.index', 'demean', 'score')
	break.up = tmp[tmp$score > thres[2],]
	break.down = tmp[tmp$score < -thres[2],]
	plot.break.up = xyplot(break.up$demean~break.up$ord.index,
						   type='p',col='red')
	plot.break.down = xyplot(break.down$demean~break.down$ord,
							 type='p',col='red')

	signal.plot = plot0 + plot.sig.up + plot.sig.down +
		plot.break.down + plot.break.up
	results = list(
		plot = signal.plot,
		break.up = break.up,
		break.down = break.down
	)
	return(signal.plot)
}



summary.zgarch = function(zgarch, trade.rule=NULL, fee.rate=0.0015, ...) {
# Summary method for strategy
# zgarch: (zgarch object), zGARCHsignal() modeling result
  df = merge(zgarch$x, zgarch$signal)
  colnames(df) = c('x', 'signal')
  if(is.null(trade.rule)) {
  	# ----- define a default trading rule -----.
  	rule.0 <- function(prev.position, new.s) {
  	# prev.position: previous position condition.
  	# new.s: new signal. (diff-cond.mean) / cond.sigma
  	  if(is.na(new.s)) {return(NA)}
	  next.pos = 0
	  if(new.s >= 2.0){
	  	next.pos = -1 # short
	  }else if(new.s <= -2.0){
	  	next.pos = 1 # long
	  }else if(abs(new.s) <= 0.6){
	  	next.pos = 0 # fill/cover
	  }else{
	  	next.pos = prev.position # do nothing
	  }
	  return(next.pos)
  	} # ---- end func ----
  	trade.rule = rule.0
  }
  df$position = Reduce(trade.rule, as.vector(df$signal), accumulate = T)
  df$position = Lag(df$position)
  df = na.omit(df)
  df$asset.ret = (df$x - Lag(df$x))
  df$realized.ret = (df$asset.ret * df$position) - fee.rate
  df = na.omit(df)
  df$equity.curve = cumsum(df$realized.ret)

  # position histograms
  df.pos = df[df$position != 0,]
  df.pos$postion = as.factor(df.pos$position)
  #hist.pos = histogram(~realized.ret|position, df.pos,
  #					 breaks=seq(from=-0.1,to=0.1,by=0.01))
  plot.signal = plotSignal.zgarch(zgarch)
  result = list(
  	df = df,
  	#hist.pos = hist.pos,
  	plot.signal = plot.signal
  )
  return(result)
}

