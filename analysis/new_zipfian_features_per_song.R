# This file expects to be run after "script5_disco_ngrams_final.R"
# and will try to extract zipfian features / measures for each individual song

disco.ngram.song <- ldply(
  sort(unique(disco.meta$song)),
  function(song){
    ## subset to this song
    ind <- which(disco.meta$song == song)
    out <- rbind(
      data.frame(type = 'pitch',
                 length = pitch.ngrams.length,
                 name = colnames(disco.pitch.dtm),
                 instances = colSums(disco.pitch.dtm[ind,]),
                 # Replaced colSums here as they fail for this class of matrix with only one row!
                 # proportion = colSums(disco.pitch.prop.dtm[ind,]) / length(ind),
                 proportion = as.numeric(disco.pitch.prop.dtm[ind,]) / length(ind),
                 songs = colSums(disco.pitch.dtm[ind,] > 0)
      ),
      data.frame(type = 'relative pitch',
                 length = pitchrel.ngrams.length,
                 name = colnames(disco.pitchrel.dtm),
                 instances = colSums(disco.pitchrel.dtm[ind,]),
                 # proportion = colSums(disco.pitchrel.prop.dtm[ind,]) / length(ind),
                 proportion = as.numeric(disco.pitchrel.prop.dtm[ind,]) / length(ind),
                 songs = colSums(disco.pitchrel.dtm[ind,] > 0)
      ),
      data.frame(type = 'rhythm',
                 length = rhythm.ngrams.length,
                 name = colnames(disco.rhythm.dtm),
                 instances = colSums(disco.rhythm.dtm[ind,]),
                 # proportion = colSums(disco.rhythm.prop.dtm[ind,]) / length(ind),
                 proportion = as.numeric(disco.rhythm.prop.dtm[ind,]) / length(ind),
                 songs = colSums(disco.rhythm.dtm[ind,] > 0)
      ),
      data.frame(type = 'relative rhythm',
                 length = rhythmrel.ngrams.length,
                 name = colnames(disco.rhythmrel.dtm),
                 instances = colSums(disco.rhythmrel.dtm[ind,]),
                 # proportion = colSums(disco.rhythmrel.prop.dtm[ind,]) / length(ind),
                 proportion = as.numeric(disco.rhythmrel.prop.dtm[ind,]) / length(ind),
                 songs = colSums(disco.rhythmrel.dtm[ind,] > 0)
      )
    )
    ## summarize frequencies of absolute/relative pitch/rhythm 2/3-grams
    out <-
      ddply(out, c('type', 'length'), function(x){
        x <- x[order(x$proportion, decreasing = TRUE),]
        x$proportion.cumulative <- cumsum(x$proportion)
        x$song.rank <- 1:nrow(x)
        return(x)
      })
    out$song <- song
    return(out)
    
})

zm.rsq.song <- list()
zm.coefs.song <- list()
## define subset of interest
for (type in c('relative pitch', 'relative rhythm')){
  for (n in 2:3){
    name <- sprintf('%s_%sgram', gsub('relative ', 'rel', type), n)
    
    N <- Ns[name]            # number of possible events
    if (is.na(N)){           #   (if unknown, estimate following zm.ll)
      N <- nrow(d)
    }
    
    d.song <- disco.ngram.song[
      disco.ngram.song$type == type &
        disco.ngram.song$length == n,
    ]
    
    ## fitted proportion predicted by zipf-mandelbrot, evaluated, but not refitted by song
    zm.rsq.song[[name]] <- list()
    zm.coefs.song[[name]] <- list()
    d.song$proportion.zm <- NA
    for (song in unique(d.song$song)){
      ind <- which(d.song$song == song)
      freqs <- d.song$instances[ind]     # counts of unique events
      class(freqs) <- 'table'  #   (let zm.ll recognize as counts)
      
      # -- 1. Use existing zipf fit across songs (this can creaty crazy R^2) --
      # normalization <- sum(1 / (1:N + zm.coefs[[name]]['b'])^zm.coefs[[name]]['s'])
      # d.song$proportion.zm[ind] <-
      #   1 / (d.song$song.rank[ind] + zm.coefs[[name]]['b'])^zm.coefs[[name]]['s'] / normalization
      # zm.rsq.song[[name]][[song]] <-
      #   1 -
      #   var(d.song$proportion[ind] - d.song$proportion.zm[ind]) /
      #   var(d.song$proportion[ind])
      
      # -- 2. Try to fit Zipf-Mandelbrot, this currently fails for quite a few cases --
      # error msg: L-BFGS-B needs finite values of 'fn'
      # tryCatch({
      #   fit <- zm.ll(freqs, N, dist = 'Zipf-Man')
      #   normalization <- sum(1 / (1:N + fit@coef['b'])^fit@coef['s'])
      #   d.song$proportion.zm[ind] <-
      #     1 / (d.song$song.rank[ind] + fit@coef['b'])^fit@coef['s'] / normalization
      #   zm.coefs.song[[name]][[song]] <- fit@coef
      #   zm.rsq.song[[name]][[song]] <-
      #     1 -
      #     var(d.song$proportion[ind] - d.song$proportion.zm[ind]) /
      #     var(d.song$proportion[ind])
      # }, error = function (cond) {
      #   warning(paste("Error occured:", name, song))
      #   warning(cond)
      #   d.song$proportion.zm[ind] <<- NA
      #   zm.coefs.song[[name]][[song]] <<- NULL
      #   zm.rsq.song[[name]][[song]] <<- NA
      # })
      
      # --  3. Fit only Zipf distribution --
      # I kept the formula here exactly the same and just replaced b with 0
      tryCatch({
        fit <- zm.ll(freqs, N, dist = 'Zipf')
        normalization <- sum(1 / (1:N + 0)^fit@coef['s'])
        d.song$proportion.zm[ind] <-
          1 / (d.song$song.rank[ind] + 0)^fit@coef['s'] / normalization
        zm.coefs.song[[name]][[song]] <- fit@coef
        zm.rsq.song[[name]][[song]] <-
          1 -
          var(d.song$proportion[ind] - d.song$proportion.zm[ind]) /
          var(d.song$proportion[ind])
      }, error = function (cond) {
        warning(paste("Error occured:", name, song))
        warning(cond)
        d.song$proportion.zm[ind] <<- NA
        zm.coefs.song[[name]][[song]] <<- NULL
        zm.rsq.song[[name]][[song]] <<- NA
      })
    }
    write.csv(d.song,
              file.path(results.dir, sprintf('disco_ngram_%s_song_final.csv', name))
    )
  }
}

zipfian_features <- data.frame(song = disco.meta$song)
for (feature in names(zm.rsq.song)) {
  zipfian_features[feature] <- as.numeric(zm.rsq.song[[feature]])
}
readr::write_csv(zipfian_features, file.path(results.dir, "disco_ngram_features.csv"))
