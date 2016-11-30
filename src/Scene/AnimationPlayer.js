const FRAMERATE = 60;
const INTERVAL_TIME = 1000 / FRAMERATE;
const debugAnimation = false;

const PLAYER_STATE = {
    STOP: 0,
    PLAY: 1,
    END: 2,
    PAUSE: 3,
};

function AnimationPlayer(dom, duration) {
    this.dom = dom;
    this.id = null;
    this.keyframe = 0;
    this.frameEvent = new CustomEvent('frameAnimation');
    this.stopEvent = new CustomEvent('stopAnimation');
    this.endEvent = new CustomEvent('endAnimation');
    this.duration = duration || 120;
    this.animation = null;
    this.resolve = null;
    this.promise = new Promise((r) => { this.resolve = r; });
    this.state = PLAYER_STATE.STOP;
}

AnimationPlayer.prototype = Object.create(AnimationPlayer.prototype);

AnimationPlayer.prototype.constructor = AnimationPlayer;

if (debugAnimation) {
    AnimationPlayer.prototype.debug = function (message) {
        /* eslint-disable no-alert, no-console */
        if (this.animation) {
            console.info('Animation ', message, ' : ', this.animation.name);
        }
        /* eslint-enable no-alert, no-console */
    };
} else {
    AnimationPlayer.prototype.debug = function () {
    };
}

AnimationPlayer.prototype.resetTimer = function () {
    if (this.id) {
        clearInterval(this.id);
        this.id = undefined;
    }

    this.keyframe = 0;
    this.duration = 0;
};

AnimationPlayer.prototype.startTimer = function () {
    this.resetTimer();
    this.id = setInterval(this.frame.bind(this), INTERVAL_TIME);
};

AnimationPlayer.prototype.start = function (duration) {
    this.debug('start');
    this.state = PLAYER_STATE.PLAY;
    this.startTimer();
    this.duration = duration === undefined ? 0 : duration;
    return new Promise((r) => { this.resolve = r; });
};

AnimationPlayer.prototype.startAnimation = function (animation) {
    this.animation = animation;
    return this.start(animation.duration);
};

AnimationPlayer.prototype.isPlayed = function () {
    return this.state === PLAYER_STATE.PLAY;
};

AnimationPlayer.prototype.isStopped = function () {
    return this.state === PLAYER_STATE.STOP;
};

AnimationPlayer.prototype.isEnded = function () {
    return this.state === PLAYER_STATE.END;
};

AnimationPlayer.prototype.finishAnimation = function () {
    this.resetTimer();
    this.animation = null;
    this.resolve();
};

AnimationPlayer.prototype.frame = function () {
    if (this.keyframe === this.duration) {
        this.debug('end');
        this.state = PLAYER_STATE.END;
        this.finishAnimation();
        this.dom.dispatchEvent(this.endEvent);
    }
    else {
        if (this.animation) {
            this.animation.animate(this.keyframe);
        }
        // console.log('frame');
        this.keyframe++;
        this.dom.dispatchEvent(this.frameEvent);
    }
};

AnimationPlayer.prototype.stop = function () {
    this.state = PLAYER_STATE.STOP;
    this.debug('stop');
    this.finishAnimation();
    this.dom.dispatchEvent(this.stopEvent);
    return this.promise;
};

AnimationPlayer.prototype.pause = function () {
    clearInterval(this.id);
    this.id = undefined;
};

AnimationPlayer.prototype.play = function () {
    this.id = setInterval(this.frame.bind(this), INTERVAL_TIME);
};


function Animation(duration, root, name) {
    this.duration = duration || 0;
    this.root = root;
    this.name = name;
}

Animation.prototype = Object.create(Animation.prototype);

Animation.prototype.animate = function (/* keyFrame*/) {

};

function AnimatedExpression(duration, root, expression, name) {
    Animation.call(this, duration, root, name);
    this.expression = expression;
}

AnimatedExpression.prototype = Object.create(Animation.prototype);

AnimatedExpression.prototype.constructor = AnimatedExpression;

AnimatedExpression.prototype.animate = function (keyFrame) {
    this.expression(this.root, keyFrame / this.duration);
};

export { Animation, AnimatedExpression };
export default AnimationPlayer;
