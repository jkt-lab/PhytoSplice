# ---- rMATS → GRanges converters ----

rmats_event_to_granges <- function(event_row, event_type) {
  
  if (event_type == "SE") {
    return(GRanges(
      seqnames = event_row$chr,
      ranges = IRanges(
        start = event_row$upstreamES,
        end   = event_row$downstreamEE
      )
    ))
  }
  
  if (event_type == "RI") {
    return(GRanges(
      seqnames = event_row$chr,
      ranges = IRanges(
        start = event_row$riExonStart_0base,
        end   = event_row$riExonEnd
      )
    ))
  }
  
  if (event_type %in% c("A3SS", "A5SS")) {
    return(GRanges(
      seqnames = event_row$chr,
      ranges = IRanges(
        start = min(event_row$longExonStart_0base,
                    event_row$shortES),
        end   = max(event_row$longExonEnd,
                    event_row$shortEE)
      )
    ))
  }
  
  if (event_type == "MXE") {
    return(GRanges(
      seqnames = event_row$chr,
      ranges = IRanges(
        start = event_row$upstreamES,
        end   = event_row$downstreamEE
      )
    ))
  }
  
  stop("Unsupported event type")
}

# ---- splicejam wrapper ----

generate_splicejam_plot <- function(
  bam_files,
  gtf,
  gr
) {
  splicejam::plotSashimi(
    bam_files = bam_files,
    gtf       = gtf,
    region    = gr
  )
}
