"use strict";
class Polynomial {
    constructor(values) {
        this.values = values;
        if (values.length === 0)
            console.log("invalid polynomial");
        this.order = values.length - 1;
    }
    dfft(coeffs, invert = false) {
        const n = coeffs.length;
        if (n === 1)
            return;
        const evenT = new Array(n / 2);
        const oddT = new Array(n / 2);
        for (let i = 0; 2 * i < n; ++i) {
            evenT[i] = coeffs[2 * i];
            oddT[i] = coeffs[2 * i + 1];
        }
        this.dfft(evenT, invert);
        this.dfft(oddT, invert);
        const pi = Math.acos(-1);
        const angle = 2 * pi / n * (invert ? -1 : 1);
        let w = { r: 1, i: 0 };
        const wn = { r: Math.cos(angle), i: Math.sin(angle) };
        for (let i = 0; 2 * i < n; ++i) {
            const t = {
                r: w.r * oddT[i].r - w.i * oddT[i].i,
                i: w.r * oddT[i].i + w.i * oddT[i].r
            };
            coeffs[i] = {
                r: evenT[i].r + t.r,
                i: evenT[i].i + t.i
            };
            coeffs[i + n / 2] = {
                r: evenT[i].r - t.r,
                i: evenT[i].i - t.i
            };
            if (invert) {
                coeffs[i].r /= 2;
                coeffs[i].i /= 2;
                coeffs[i + n / 2].r /= 2;
                coeffs[i + n / 2].i /= 2;
            }
            w = {
                r: w.r * wn.r - w.i * wn.i,
                i: w.r * wn.i + w.i * wn.r
            };
        }
    }
    prune(array) {
        let i = 0;
        while (array.length > 1) {
            if (array[i] !== 0)
                break;
            array.splice(i, 1);
        }
    }
    multiply(poly) {
        let n = 1;
        while (n < poly.values.length + poly.values.length)
            n <<= 1;
        const fftA = new Array(n);
        const fftB = new Array(n);
        for (let i = 0; i < n; ++i) {
            fftA[i] = { r: i < this.values.length ? this.values[i] : 0, i: 0 };
            fftB[i] = { r: i < poly.values.length ? poly.values[i] : 0, i: 0 };
        }
        this.dfft(fftA);
        this.dfft(fftB);
        for (let i = 0; i < n; ++i) {
            fftA[i] = {
                r: fftA[i].r * fftB[i].r - fftA[i].i * fftB[i].i,
                i: fftA[i].r * fftB[i].i + fftA[i].i * fftB[i].r
            };
        }
        this.dfft(fftA, true);
        const result = new Array(n);
        for (let i = 0; i < n; ++i) {
            result[i] = Math.round(fftA[i].r);
        }
        while (result.length >= (this.values.length + poly.values.length)) {
            result.pop();
        }
        return new Polynomial(result);
    }
    add(poly) {
        let lmax = Math.max(poly.order, this.order);
        let accumulator = new Array(lmax + 1).fill(0);
        let offset1 = 0, offset2 = 0;
        for (let i = 0; i <= lmax; i++) {
            if (lmax <= this.order + i)
                accumulator[i] += this.values[i - offset1];
            else
                offset1 += 1;
            if (lmax <= poly.order + i)
                accumulator[i] += poly.values[i - offset2];
            else
                offset2 += 1;
        }
        this.prune(accumulator);
        return new Polynomial(accumulator);
    }
    sub(poly) {
        let lmax = Math.max(poly.order, this.order);
        let accumulator = new Array(lmax + 1).fill(0);
        let offset1 = 0, offset2 = 0;
        for (let i = 0; i <= lmax; i++) {
            if (lmax <= this.order + i)
                accumulator[i] += this.values[i - offset1];
            else
                offset1 += 1;
            if (lmax <= poly.order + i)
                accumulator[i] -= poly.values[i - offset2];
            else
                offset2 += 1;
        }
        this.prune(accumulator);
        return new Polynomial(accumulator);
    }
    calculate(x) {
        if (this.order === 0)
            return this.values[0];
        let accumulator = 0;
        for (let i = 0; i < this.values.length; i++) {
            accumulator += this.values[i] * (x ** (this.order - i));
        }
        return accumulator;
    }
    derive() {
        if (this.order === 0)
            return new Polynomial([0]);
        let vs = [];
        for (let i = 0; i < this.values.length - 1; i++) {
            vs.push(this.values[i] * (this.order - i));
        }
        return new Polynomial(vs);
    }
    estimateRoot(margin = 0.005) {
        if (this.order === 1)
            return { x: -this.values[1] / this.values[0], value: 0 };
        if (this.order === 0)
            return null;
        let firstDer = this.derive();
        let secondDer = firstDer.derive();
        let x = 0;
        let loss = 10 ** 10, bestLoss = 10 ** 10, bestX = 0;
        let accumulator = 0;
        let stop = false;
        let count = 1;
        while (!stop) {
            console.log(count++);
            let bv = this.calculate(x), fv = firstDer.calculate(x), sv = secondDer.calculate(x);
            x -= bv * fv / (fv ** 2 - bv * sv / 2);
            loss = this.calculate(x);
            if (Math.abs(loss) < Math.abs(bestLoss)) {
                bestLoss = loss;
                bestX = x;
            }
            else {
                accumulator += 1;
            }
            if (accumulator > 5 || Math.abs(loss) < margin) {
                stop = true;
                break;
            }
        }
        return {
            x: +bestX.toFixed(5),
            value: +loss.toFixed(5)
        };
    }
}
