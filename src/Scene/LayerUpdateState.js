const UPDATE_IDLE = 0;
const UPDATE_PENDING = 1;
const UPDATE_ERROR = 2;

const PAUSE_BETWEEN_ERRORS = [ 1000.0, 3000.0, 7000.0, 60000.0 ];

function LayerUpdateState() {
    this.state = UPDATE_IDLE;
    this.lastErrorTimestamp = 0;
    this.errorCount = 0;

}

LayerUpdateState.prototype.constructor = LayerUpdateState;

LayerUpdateState.prototype.canTryUpdate = function(timestamp) {
    switch (this.state) {
        case UPDATE_IDLE: {
            return true;
        }
        case UPDATE_PENDING: {
            return false;
        }
        case UPDATE_ERROR:
        default: {
            let errorDuration = this.waitDurationUntilNextTry();
            return errorDuration <= (timestamp - this.lastErrorTimestamp);
        }
    }
}

LayerUpdateState.prototype.waitDurationUntilNextTry = function() {
    if (this.state != UPDATE_ERROR) {
        return 0;
    }
    let idx =
        Math.max(0, Math.min(this.errorCount, PAUSE_BETWEEN_ERRORS.length) - 1);

    return PAUSE_BETWEEN_ERRORS[idx];
}

LayerUpdateState.prototype.try = function() {
    this.state = UPDATE_PENDING;
}

LayerUpdateState.prototype.success = function() {
    this.lastErrorTimestamp = 0;
    this.state = UPDATE_IDLE;
}

LayerUpdateState.prototype.failure = function(timestamp) {
    this.lastErrorTimestamp = timestamp;
    this.state = UPDATE_ERROR;
    this.errorCount++;
}

export default LayerUpdateState;
