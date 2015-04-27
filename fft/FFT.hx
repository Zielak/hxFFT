package fft;

/**
 * Performs an in-place complex FFT.
 *
 * Released under the MIT License
 *
 * Copyright (c) 2010 Gerald T. Beauregard
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

import haxe.ds.Vector;

class FFT
{
    public static var FORWARD:Bool = false;
    public static var INVERSE:Bool = true;

    private var m_logN:Int = 0;            // log2 of FFT size
    private var m_N:Int = 0;               // FFT size
    private var m_invN:Float;              // Inverse of FFT length

    private var m_X:Vector<FFTElement>;    // Vector of linked list elements

    /**
     *
     */
    public function new()
    {
    }

    /**
     * Initialize class to perform FFT of specified size.
     *
     * @param   logN    Log2 of FFT length. e.g. for 512 pt FFT, logN = 9.
     */
    public function init( logN:Int ):Void
    {
        m_logN = logN;
        m_N = 1 << m_logN;
        m_invN = 1.0/m_N;

        // Allocate elements for linked list of complex numbers.
        m_X = new Vector<FFTElement>(m_N);
        for( k in 0...m_N ){
            m_X[k] = new FFTElement();
        }

        // Set up "next" pointers.
        var _mNone = m_N-1;
        for( k in 0..._mNone ){
            m_X[k].next = m_X[k+1];
        }

        // Specify target for bit reversal re-ordering.
        for( k in 0...m_N ){
            m_X[k].revTgt = BitReverse(k,logN);
        }
    }

    /**
     * Performs in-place complex FFT.
     *
     * @param   xRe     Real part of input/output
     * @param   xIm     Imaginary part of input/output
     * @param   inverse If true (INVERSE), do an inverse FFT
     */
    public function run(
        xRe:Vector<Float>,
        xIm:Vector<Float>,
        inverse:Bool = false ):Void
    {
        var numFlies:Int = m_N >> 1;   // Number of butterflies per sub-FFT
        var span:Int = m_N >> 1;       // Width of the butterfly
        var spacing:Int = m_N;         // Distance between start of sub-FFTs
        var wIndexStep:Int = 1;        // Increment for twiddle table index

        // Copy data into linked complex number objects
        // If it's an IFFT, we divide by N while we're at it
        var x:FFTElement = m_X[0];
        var k:Int = 0;
        var scale:Float = inverse ? m_invN : 1.0;

        while (x)
        {
            x.re = scale*xRe[k];
            x.im = scale*xIm[k];
            x = x.next;
            k++;
        }

        // For each stage of the FFT
        for( stage in 0...m_logN )
        {
            // Compute a multiplier factor for the "twiddle factors".
            // The twiddle factors are complex unit vectors spaced at
            // regular angular intervals. The angle by which the twiddle
            // factor advances depends on the FFT stage. In many FFT
            // implementations the twiddle factors are cached, but because
            // vector lookup is relatively slow in ActionScript, it's just
            // as fast to compute them on the fly.
            var wAngleInc:Float = wIndexStep * 2.0*Math.PI/m_N;
            if ( inverse == false )
            {
                // Corrected 3 Aug 2011. Had this condition backwards before, so FFT was IFFT, and vice-versa!
                wAngleInc *= -1;
            }

            var wMulRe:Float = Math.cos(wAngleInc);
            var wMulIm:Float = Math.sin(wAngleInc);

            var start:Int = 0;
            while( start < m_N )
            {
                start += spacing;

                var xTop:FFTElement = m_X[start];
                var xBot:FFTElement = m_X[start+span];

                var wRe:Float = 1.0;
                var wIm:Float = 0.0;

                // For each butterfly in this stage
                for( flyCount in 0...numFlies )
                {
                    // Get the top & bottom values
                    var xTopRe:Float = xTop.re;
                    var xTopIm:Float = xTop.im;
                    var xBotRe:Float = xBot.re;
                    var xBotIm:Float = xBot.im;

                    // Top branch of butterfly has addition
                    xTop.re = xTopRe + xBotRe;
                    xTop.im = xTopIm + xBotIm;

                    // Bottom branch of butterly has subtraction,
                    // followed by multiplication by twiddle factor
                    xBotRe = xTopRe - xBotRe;
                    xBotIm = xTopIm - xBotIm;
                    xBot.re = xBotRe*wRe - xBotIm*wIm;
                    xBot.im = xBotRe*wIm + xBotIm*wRe;

                    // Advance butterfly to next top & bottom positions
                    xTop = xTop.next;
                    xBot = xBot.next;

                    // Update the twiddle factor, via complex multiply
                    // by unit vector with the appropriate angle
                    // (wRe + j wIm) = (wRe + j wIm) x (wMulRe + j wMulIm)
                    var tRe:Float = wRe;
                    wRe = wRe*wMulRe - wIm*wMulIm;
                    wIm = tRe*wMulIm + wIm*wMulRe;
                }
            }

            numFlies >>= 1;   // Divide by 2 by right shift
            span >>= 1;
            spacing >>= 1;
            wIndexStep <<= 1;     // Multiply by 2 by left shift
        }

        // The algorithm leaves the result in a scrambled order.
        // Unscramble while copying values from the complex
        // linked list elements back to the input/output vectors.
        x = m_X[0];
        while (x)
        {
            var target:Int = x.revTgt;
            xRe[target] = x.re;
            xIm[target] = x.im;
            x = x.next;
        }
    }

    /**
     * Do bit reversal of specified number of places of an int
     * For example, 1101 bit-reversed is 1011
     *
     * @param   x       Number to be bit-reverse.
     * @param   numBits Number of bits in the number.
     */
    private function BitReverse(
        x:Int,
        numBits:Int):Int
    {
        var y:Int = 0;
        for( i in 0...numBits )
        {
            y <<= 1;
            y |= x & 0x0001;
            x >>= 1;
        }
        return y;
    }
}

class FFTElement
{
    public var re:Float;          // Real component
    public var im:Float;          // Imaginary component
    public var next:FFTElement;   // Next element in linked list
    public var revTgt:Int;        // Target position post bit-reversal

    public function new(){
        re = 0;
        im = 0;
        next = null;
        revTgt = 0;
    }
}
