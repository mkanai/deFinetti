#' Plot a de Finetti diagram
#' 
#' @param x a numeric \code{matrix}
#' @param file a string of output filename
#' @param main a string of axis label
#' @param vertexlab a string \code{vector} of vertex labels
#' @param cex a numeric of character expansion factor
#' @param cex.main a numeric of magnification for titles relative to \code{cex}
#' @param cex.lab a numeric of magnification for labels relative to \code{cex}
#' @param pch a symbol for plotting markers
#' @param markerlab a string \code{vector} of marker labels
#' @param markerpos a position specifier for marker labels
#' @param markercex a numeric of magnification for markers relative to \code{cex}
#' @param markercol a color for plotting marker symbols
#' @param markerbg a background (fill) color for plotting marker symbols given by \code{pch = 21:25}
#' @param hwcurve a logical indicating whether to draw a curve of Hardy-Weinberg Equilibrium
#' @param with_F_color a logical indicating whether to draw colors based on F-statistics
#' @param without_markers a logical indicating whether to plot markers
#' @param pdf_out a logical indicating whether to output a diagram
#' @export
deFinetti.plot <- function(x, file = "", main = "", vertexlab = colnames(x), 
                           cex = 1, cex.main = 2, cex.lab = 2, pch = 19,                           
                           markerlab = NULL, markerpos = 1, markercex = 1, markercol = "black", markerbg = "black", 
                           hwcurve = TRUE, with_F_color = FALSE, without_markers = FALSE, pdf_out = FALSE, ...) {
    # Graphics
    oldpar = par(no.readonly = T)
    on.exit(par(oldpar))
    par(pty = "m")
    par(xpd = TRUE)
    if (main != "") {
        par(oma = c(0, 0, 4, 0))
    }    

    # Init
    if (pdf_out) {
        pdf(ifelse(file != "", file, "Definetti_Diagram.pdf"), paper = "a4r")
    }
    
    if (without_markers) {
        aa <- c(100, 100, 100)
        x <- data.frame(AA = aa, AB = aa, BB = aa)
    }
    
    if (is.vector(x)) {
        if (length(x) != 3) {
            stop("x should have three columns.")
        }
        x <- matrix(x, ncol = 3, dimnames = list(c("1", names(x))))
    } else {
        x <- as.matrix(x)
    }
    
    if (ncol(x) != 3) {
        stop("x should have three columns")
    }
    
    if (any(x < 0)) {
        stop("x should be non-negative value")
    }
    
    if (sum(apply(x, 1, sum)) == nrow(x)) {
        xcom <- x
    } else {
        xr <- x
        if (nrow(x) == 1) {
            xcom <- x/sum(x)
        } else {
            dr <- diag(apply(x, 1, sum))
            xcom <- solve(dr) %*% x
        }
    }
    
    # Legend
    if (with_F_color) {        
        layout(matrix(c(1, rep(2, 3)), 1, 4))
        plot.new()        
        xpos <- .5
        ypos <- 0
        width <- .1
        cols <- colorRamp(c("#ff0000", "white", "#0000ff"))
        rect(xpos, ypos+(0:999)/1000, xpos + width, ypos + (1:1000)/1000, col = rgb(cols(0:999/999)/255), border = NA)
        text(xpos, 0, expression(italic(F) == 1), pos = 2, cex = cex*cex.lab)
        text(xpos, .5, expression(italic(F) == 0), pos = 2, cex = cex*cex.lab)
        text(xpos, 1, expression(italic(F) == -1), pos = 2, cex = cex*cex.lab)
    }
        
    # Triangle
    s <- 1/sqrt(3)
    M <- matrix(c(-s, 0, 0, 1, s, 0), ncol = 2, byrow = T)
    M2 <- matrix(c(-s, 0, 0, 1, -2*s, 1, -2*s, 0), ncol = 2, byrow = T)
    M3 <- matrix(c(s, 0, 0, 1, 3*s, 1, 3*s, 0), ncol = 2, byrow = T)
        
    plot(M[,1], M[,2], type = "n", axes = F, xlab = "", ylab = "", pch = 19, asp = 1, ...)    
    polygon(M, col = NULL)
    
    # F Color
    if (with_F_color) {
        for (step in 0:255) {
            color_step <- 255 - step
            f <- color_step / 255
            
            p <- seq(0, 1, by = .005)
            HWb <- cbind(f*p + (p^2)*(1-f), 2*(1-f)*p*(1-p), f*(1-p) + (1-f)*(1-p)^2) %*% M            
            HWr <- cbind(-f*p + (p^2)*(1+f), 2*(1+f)*p*(1-p), -f*(1-p) + (1+f)*(1-p)^2) %*% M
            
            points(HWb[,1], HWb[,2], type = "l", col = rgb(1, 1-f, 1-f), lwd = 4)
            points(HWr[,1], HWr[,2], type = "l", col = rgb(1-f, 1-f, 1), lwd = 4)
        }
        
        # Over-draw
        polygon(M2, col = "white", border = "white", lwd = 4)
        polygon(M3, col = "white", border = "white", lwd = 4)
        
        lines(c(-s, 0, s, -s), c(0, 1, 0, 0), type = "l", col = "black", lwd = 4)
    }
    
    # HW Curve
    if (hwcurve) {
        p <- seq(0, 1, by = .005)
        HWc <- cbind(p^2, 2*p*(1-p), (1-p)^2) %*% M
        points(HWc[,1], HWc[,2], type = "l", col = "black")
    }

    # Markers
    if(!without_markers) {
        xc <- xcom %*% M
        points(xc[, 1], xc[, 2], pch = pch, bg = markerbg, col = markercol, cex = markercex)
        if (!is.null(markerlab)) {text(xc[, 1], xc[, 2], markerlab, cex = cex*cex.lab, pos = markerpos)}
    }

    # Vertex Labels
    eps <- .03 * cex*cex.lab
    Mlab <- M + matrix(c(-eps, 0, 0, eps, eps, 0), ncol = 2, byrow = T)
    text(Mlab[,1], Mlab[,2], vertexlab, cex = cex*cex.lab)
    
    # Title
    if (main != "") {
        mtext(side = 3, line = 1, outer = T, text = main, cex = cex*cex.main)
    }
    
    if (pdf_out) {dev.off()}
}